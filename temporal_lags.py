#!/usr/bin/env python3
"""
Temporal separation between the pre-event and post-event acquisitions,
per burn event and for the oil site, and whether separability depends
on it.

Half the features in the analysis are differences between a pre-event
and a post-event scene. If the fires are observed across much longer
intervals than the spill, part of what the differenced features measure
could be phenology rather than burning. This quantifies the exposure.

Reported per event:
    dt_pre   days from the pre-event scene to the start of the event
    dt_post  days from the end of the event to the post-event scene
    dt_tot   days between the two scenes

and then the rank correlation of each against the event's own
separability from oil on the features the manuscript discusses. A
correlation near zero says the differenced features are not being
driven by how long the interval happened to be.

Dependencies: numpy only.

Usage:
    python3 temporal_lags.py
"""

import csv
import datetime as dt
import os
import numpy as np

from oil_burn_separability import (
    load, subset, column, separability, _midranks, MIN_EVENT_N,
    OUT_RESULTS, pretty,
)

REPORT = ['d_NBR', 'd_FMR', 'd_logBAI', 'd_BAI', 'MIRBI', 'NBR2']


def parse(s):
    s = (s or '').strip()[:10]
    for fmt in ('%Y-%m-%d', '%Y/%m/%d', '%d/%m/%Y'):
        try:
            return dt.datetime.strptime(s, fmt).date()
        except ValueError:
            pass
    return None


def spearman(x, y):
    """Rank correlation, mid-ranks for ties. Pearson on the ranks."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan
    rx, ry = _midranks(x[m]), _midranks(y[m])
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else np.nan


def describe(name, a):
    a = a[np.isfinite(a)]
    if not a.size:
        return
    print('  %-8s n=%2d  median %5.1f  IQR %5.1f-%5.1f  range %5.1f-%5.1f'
          % (name, a.size, np.median(a), np.percentile(a, 25),
             np.percentile(a, 75), a.min(), a.max()))


def main():
    path = 'oil_burn_samples.csv'
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)
    rows, feats = load(path)
    oil = subset(rows, 'oil')
    burn = subset(rows, 'burn')

    print('=' * 70)
    print('TEMPORAL SEPARATION OF THE IMAGE PAIRS')
    print('=' * 70)

    # --- the oil site -------------------------------------------------
    r = oil[0]
    o_pre, o_post = parse(r['pre_date']), parse(r['post_date'])
    o_ev, o_end = parse(r['event_date']), parse(r['event_end_date'])
    if o_ev is None:
        o_ev = o_end
    if o_end is None:
        o_end = o_ev
    print('oil  pre %s  event %s..%s  post %s'
          % (r['pre_date'][:10], r['event_date'][:10],
             r['event_end_date'][:10], r['post_date'][:10]))
    if None not in (o_pre, o_post, o_ev, o_end):
        print('     dt_pre %d d   dt_post %d d   dt_tot %d d'
              % ((o_ev - o_pre).days, (o_post - o_end).days,
                 (o_post - o_pre).days))
    print()

    # --- burn events ---------------------------------------------------
    big, meta = {}, {}
    for e in sorted({x['event_id'] for x in burn}):
        b = [x for x in burn if x['event_id'] == e]
        if len(b) >= MIN_EVENT_N:
            big[e] = b
            meta[e] = b[0]

    recs = []
    for e, b in big.items():
        m = meta[e]
        pre, post = parse(m['pre_date']), parse(m['post_date'])
        ev, end = parse(m['event_date']), parse(m['event_end_date'])
        if ev is None:
            ev = end
        if end is None:
            end = ev
        rec = {
            'event_id': e,
            'n_pixels': len(b),
            'area_ha': m.get('area_ha', ''),
            'pre_date': m['pre_date'][:10],
            'post_date': m['post_date'][:10],
            'dt_pre': (ev - pre).days if None not in (ev, pre) else np.nan,
            'dt_post': (post - end).days if None not in (post, end) else np.nan,
            'dt_tot': (post - pre).days if None not in (post, pre) else np.nan,
            'dt_event': (end - ev).days if None not in (end, ev) else np.nan,
        }
        recs.append(rec)

    for k in ('dt_pre', 'dt_post', 'dt_tot', 'dt_event'):
        describe(k, np.array([r[k] for r in recs], dtype=float))
    print()

    # --- separability per event, per feature ---------------------------
    oil_c = {f: column(oil, f) for f in feats}
    for f in REPORT:
        if f not in feats:
            continue
        for rec in recs:
            b = column(big[rec['event_id']], f)
            rec['S_%s' % f] = separability(oil_c[f], b)

    print('RANK CORRELATION OF EVENT SEPARABILITY WITH THE INTERVAL')
    print('%-12s %8s %8s %8s' % ('feature', 'dt_pre', 'dt_post', 'dt_tot'))
    print('-' * 40)
    for f in REPORT:
        k = 'S_%s' % f
        if k not in recs[0]:
            continue
        s = np.array([r[k] for r in recs], dtype=float)
        print('%-12s %8.3f %8.3f %8.3f' % (
            pretty(f),
            spearman(np.array([r['dt_pre'] for r in recs], dtype=float), s),
            spearman(np.array([r['dt_post'] for r in recs], dtype=float), s),
            spearman(np.array([r['dt_tot'] for r in recs], dtype=float), s)))
    print()

    os.makedirs(OUT_RESULTS, exist_ok=True)
    p = os.path.join(OUT_RESULTS, 'event_intervals.csv')
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)
    print('written: %s  (%d events)' % (p, len(recs)))


if __name__ == '__main__':
    main()
