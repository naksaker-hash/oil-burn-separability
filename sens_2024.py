#!/usr/bin/env python3
"""
2024 sensitivity run.

18 of the 58 matched burn events fall in 2024, contributing 4167 of the
14726 burn pixels. If the dissociation is an artefact of that one fire
season, dropping those events should collapse it. This script runs the
headline analysis twice -- all events, then 2024 excluded -- and prints
the two side by side.

Every statistic is imported from oil_burn_separability.py, unchanged, so
the two runs cannot drift apart by reimplementation. Only the row set
differs.

The oil class is untouched: the plume is 2018 and there is only one.
Excluding 2024 removes burn events and their own matched control rings.

Usage:
    python3 sens_2024.py [path/to/oil_burn_samples.csv]
"""

import csv
import os
import sys
import numpy as np

import oil_burn_separability as ob

OUT_RESULTS = 'results'

# The features named in the manuscript's dissociation table, in the
# order they appear there. Kept explicit so the comparison table has a
# stable shape whatever the re-ranking does.
HEADLINE = ['d_NBR', 'd_FMR', 'NBR', 'd_logBAI', 'd_BAI',
            'MIRBI', 'NBR2', 'NDVI', 'B11', 'ARVI']

# Event-cluster bootstrap is the expensive step, so it is restricted to
# the features the argument actually turns on.
CI_FEATURES = ['MIRBI', 'NBR2', 'NDVI', 'd_NBR', 'd_FMR', 'd_logBAI']

EXCLUDE_YEAR = '2024'


def year(row):
    return row['event_date'][:4]


def analyse(oil, edge, clean, burn, ctrl, feats, label):
    """The headline numbers for one row set. Mirrors main() in the
    parent script, minus the figures and the printing."""
    res = {'label': label}

    ev = sorted({r['event_id'] for r in burn})
    res['n_events'] = len(ev)
    res['n_burn_px'] = len(burn)

    oil_c = {f: ob.column(oil, f) for f in feats}
    clean_c = {f: ob.column(clean, f) for f in feats}
    burn_c = {f: ob.column(burn, f) for f in feats}

    big_events = {}
    for e in ev:
        b = [r for r in burn if r['event_id'] == e]
        if len(b) >= ob.MIN_EVENT_N:
            big_events[e] = b
    res['n_big_events'] = len(big_events)
    ev_cols = {f: {k: ob.column(v, f) for k, v in big_events.items()}
               for f in feats}

    # ---- label quality, burn vs its own control ring on d_NBR ----
    qc = []
    for e in ev:
        b = [r for r in burn if r['event_id'] == e]
        c = [r for r in ctrl if r['event_id'] == e]
        a = ob.auc(ob.column(c, 'd_NBR'), ob.column(b, 'd_NBR'))
        if np.isfinite(a):
            qc.append(a)
    qc = np.array(qc)
    res['qc_n'] = int(qc.size)
    res['qc_median'] = float(np.median(qc)) if qc.size else np.nan
    res['qc_ge_090'] = int((qc >= 0.90).sum()) if qc.size else 0

    # ---- per-feature separability ----
    table = []
    for f in feats:
        o, c, b = oil_c[f], clean_c[f], burn_c[f]
        s_clean = ob.separability(o, c)
        s_burn = ob.separability(o, b)
        if not (np.isfinite(s_clean) and np.isfinite(s_burn)):
            continue
        per_ev = np.array([s for s in
                           (ob.separability(o, a) for a in ev_cols[f].values())
                           if np.isfinite(s)])
        table.append({
            'feature': f,
            'sep_oil_clean': s_clean,
            'sep_oil_burn': s_burn,
            'per_event_median': np.median(per_ev) if per_ev.size else np.nan,
            'per_event_n_gt_080': int((per_ev > 0.80).sum()) if per_ev.size else 0,
            'per_event_n_lt_030': int((per_ev < 0.30).sum()) if per_ev.size else 0,
            'per_event_n': int(per_ev.size),
        })
    table.sort(key=lambda d: -d['sep_oil_clean'])
    by_burn = sorted(table, key=lambda d: -d['sep_oil_burn'])
    res['table'] = {d['feature']: d for d in table}
    res['rank_clean'] = [d['feature'] for d in table]
    res['rank_burn'] = [d['feature'] for d in by_burn]

    sc = np.array([d['sep_oil_clean'] for d in table])
    sb = np.array([d['sep_oil_burn'] for d in table])
    res['n_feats'] = len(table)
    res['median_clean'] = float(np.median(sc))
    res['median_burn'] = float(np.median(sb))
    res['n_gt090_clean'] = int((sc > 0.90).sum())
    res['n_gt090_burn'] = int((sb > 0.90).sum())

    # ---- event-cluster bootstrap for the features that matter ----
    res['ci'] = {}
    for f in CI_FEATURES:
        if f in ev_cols:
            res['ci'][f] = ob.boot_ci_events(oil_c[f], ev_cols[f])

    # ---- substrate-matched sensitivity ----
    plume_ndvi = np.nanmedian(ob.column(oil, 'pre_NDVI_px'))
    plume_slope = np.nanmedian(ob.column(oil, 'slope_px'))
    bn = ob.column(burn, 'pre_NDVI_px')
    bs = ob.column(burn, 'slope_px')
    keep = (np.abs(bn - plume_ndvi) <= ob.NDVI_TOL) & \
           (np.abs(bs - plume_slope) <= ob.SLOPE_TOL)
    burn_m = [r for r, k in zip(burn, keep) if k]
    res['n_matched'] = len(burn_m)
    res['matched'] = {}
    if len(burn_m) >= ob.MIN_N:
        for f in HEADLINE:
            if f in oil_c:
                res['matched'][f] = ob.separability(oil_c[f],
                                                    ob.column(burn_m, f))

    # ---- two-feature discriminant, leave-one-burn-event-out ----
    det = table[0]['feature']
    dis = by_burn[0]['feature']
    res['det'], res['dis'] = det, dis
    rng = np.random.default_rng(ob.BOOT_SEED)
    o2 = np.column_stack([oil_c[det], oil_c[dis]])
    o2 = o2[ob.clean_rows([o2[:, 0], o2[:, 1]])]
    ev2 = {}
    for k in big_events:
        m = np.column_stack([ev_cols[det][k], ev_cols[dis][k]])
        m = m[ob.clean_rows([m[:, 0], m[:, 1]])]
        if len(m) >= ob.MIN_N:
            ev2[k] = m
    keys = list(ev2)
    scores = []
    for _ in range(ob.LOEO_REPEATS):
        idx = rng.permutation(len(o2))
        half = len(o2) // 2
        otr, ote = o2[idx[:half]], o2[idx[half:]]
        for k in keys:
            tr = [ev2[kk] for kk in keys if kk != k]
            if not tr:
                continue
            btr = np.vstack(tr)
            bte = ev2[k]
            if len(btr) < ob.MIN_N or len(bte) < ob.MIN_N:
                continue
            w = ob.fisher_ld(btr, otr)
            s = ob.separability(ote @ w, bte @ w)
            if np.isfinite(s):
                scores.append(s)
    scores = np.array(scores)
    res['loeo_n'] = int(scores.size)
    if scores.size:
        res['loeo_median'] = float(np.median(scores))
        res['loeo_p10'] = float(np.percentile(scores, 10))
        res['loeo_p90'] = float(np.percentile(scores, 90))
    else:
        res['loeo_median'] = res['loeo_p10'] = res['loeo_p90'] = np.nan
    return res


def fmt(v, nd=3):
    return 'n/a' if v is None or not np.isfinite(v) else '%.*f' % (nd, v)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else ob.CSV_DEFAULT
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)
    os.makedirs(OUT_RESULTS, exist_ok=True)

    rows, feats = ob.load(path)
    oil = ob.subset(rows, 'oil')
    edge = ob.subset(rows, 'oil_edge')
    clean = ob.subset(rows, 'control_clean')
    burn = ob.subset(rows, 'burn')
    ctrl = ob.subset(rows, 'control_unburned')

    burn_x = [r for r in burn if year(r) != EXCLUDE_YEAR]
    ctrl_x = [r for r in ctrl if year(r) != EXCLUDE_YEAR]

    print('=' * 70)
    print('2024 SENSITIVITY RUN')
    print('=' * 70)
    print('source : %s' % path)
    print('dropped: %d burn px, %d control px, %d events'
          % (len(burn) - len(burn_x), len(ctrl) - len(ctrl_x),
             len({r['event_id'] for r in burn})
             - len({r['event_id'] for r in burn_x})))
    print()

    full = analyse(oil, edge, clean, burn, ctrl, feats, 'all events')
    print('run 1 of 2 done (all events)')
    excl = analyse(oil, edge, clean, burn_x, ctrl_x, feats, 'no 2024')
    print('run 2 of 2 done (2024 excluded)')
    print()

    print('SAMPLE')
    print('%-28s %12s %12s' % ('', 'all events', 'no 2024'))
    print('-' * 54)
    for k, lab in [('n_events', 'burn events'),
                   ('n_burn_px', 'burn pixels'),
                   ('n_big_events', 'events >= 100 px'),
                   ('qc_n', 'events scored for labels')]:
        print('%-28s %12d %12d' % (lab, full[k], excl[k]))
    print('%-28s %12s %12s' % ('label quality, median AUC',
                               fmt(full['qc_median']), fmt(excl['qc_median'])))
    print('%-28s %12d %12d' % ('  of those >= 0.90',
                               full['qc_ge_090'], excl['qc_ge_090']))
    print()

    print('SEPARABILITY ACROSS ALL %d FEATURES' % full['n_feats'])
    print('%-28s %12s %12s' % ('', 'all events', 'no 2024'))
    print('-' * 54)
    print('%-28s %12s %12s' % ('median vs clean soil',
                               fmt(full['median_clean']), fmt(excl['median_clean'])))
    print('%-28s %12s %12s' % ('median vs burn',
                               fmt(full['median_burn']), fmt(excl['median_burn'])))
    print('%-28s %12d %12d' % ('> 0.90 vs clean soil',
                               full['n_gt090_clean'], excl['n_gt090_clean']))
    print('%-28s %12d %12d' % ('> 0.90 vs burn',
                               full['n_gt090_burn'], excl['n_gt090_burn']))
    print()

    print('HEADLINE FEATURES, separability vs burn')
    print('%-10s %9s %9s %8s   %-22s' %
          ('feature', 'all', 'no 2024', 'change', '95% CI, no 2024'))
    print('-' * 64)
    for f in HEADLINE:
        a = full['table'].get(f)
        b = excl['table'].get(f)
        if not a or not b:
            continue
        ci = excl['ci'].get(f)
        cis = ('[%s, %s]' % (fmt(ci[0]), fmt(ci[1]))) if ci else ''
        print('%-10s %9.3f %9.3f %+8.3f   %-22s' %
              (f, a['sep_oil_burn'], b['sep_oil_burn'],
               b['sep_oil_burn'] - a['sep_oil_burn'], cis))
    print()

    print('PER BURN EVENT (events with >= 100 burn pixels)')
    print('%-10s %11s %11s %11s %11s' %
          ('feature', 'median all', 'median no24', '>0.80 all', '>0.80 no24'))
    print('-' * 58)
    for f in ['MIRBI', 'NBR2', 'NDVI', 'd_NBR', 'd_logBAI']:
        a = full['table'].get(f)
        b = excl['table'].get(f)
        if not a or not b:
            continue
        print('%-10s %11s %11s %6d/%-4d %6d/%-4d' %
              (f, fmt(a['per_event_median']), fmt(b['per_event_median']),
               a['per_event_n_gt_080'], a['per_event_n'],
               b['per_event_n_gt_080'], b['per_event_n']))
    dl_a = full['table'].get('d_logBAI')
    dl_b = excl['table'].get('d_logBAI')
    if dl_a and dl_b:
        print('d_logBAI < 0.30 on %d of %d events (all), %d of %d (no 2024)'
              % (dl_a['per_event_n_lt_030'], dl_a['per_event_n'],
                 dl_b['per_event_n_lt_030'], dl_b['per_event_n']))
    print()

    print('SUBSTRATE-MATCHED BURN PIXELS')
    print('  matched pixels: %d (all), %d (no 2024)'
          % (full['n_matched'], excl['n_matched']))
    print('%-10s %11s %11s %11s %11s' %
          ('feature', 'all burn', 'all match', 'no24 burn', 'no24 match'))
    print('-' * 58)
    for f in ['MIRBI', 'NBR2', 'NDVI']:
        if f not in full['matched'] or f not in excl['matched']:
            continue
        print('%-10s %11s %11s %11s %11s' %
              (f, fmt(full['table'][f]['sep_oil_burn']), fmt(full['matched'][f]),
               fmt(excl['table'][f]['sep_oil_burn']), fmt(excl['matched'][f])))
    print()

    print('TWO-FEATURE DISCRIMINANT, leave-one-burn-event-out')
    print('%-28s %12s %12s' % ('', 'all events', 'no 2024'))
    print('-' * 54)
    print('%-28s %12s %12s' % ('pair',
                               '%s+%s' % (full['det'], full['dis']),
                               '%s+%s' % (excl['det'], excl['dis'])))
    print('%-28s %12d %12d' % ('folds', full['loeo_n'], excl['loeo_n']))
    print('%-28s %12s %12s' % ('median separability',
                               fmt(full['loeo_median']), fmt(excl['loeo_median'])))
    print('%-28s %12s %12s' % ('10th percentile',
                               fmt(full['loeo_p10']), fmt(excl['loeo_p10'])))
    print('%-28s %12s %12s' % ('90th percentile',
                               fmt(full['loeo_p90']), fmt(excl['loeo_p90'])))
    print()

    # ---- year jackknife on the discriminating index ----------------
    # 2024 is the year the reviewer will ask about because it is the
    # largest, but the question it stands for is whether ANY single
    # season carries the result. Dropping each year in turn answers that
    # at almost no extra cost.
    print('YEAR JACKKNIFE, %s separability vs burn' % excl['dis'])
    print('%-10s %9s %9s %9s' % ('dropped', 'events', 'burn px', 'sep'))
    print('-' * 40)
    f = excl['dis']
    o = ob.column(oil, f)
    years = sorted({year(r) for r in burn})
    jack = []
    for y in years:
        keep = [r for r in burn if year(r) != y]
        nev = len({r['event_id'] for r in keep})
        s = ob.separability(o, ob.column(keep, f))
        jack.append((y, s))
        print('%-10s %9d %9d %9s' % (y, nev, len(keep), fmt(s)))
    js = np.array([s for _, s in jack if np.isfinite(s)])
    if js.size:
        print('range over jackknife: %.3f to %.3f (full sample %.3f)'
              % (js.min(), js.max(), full['table'][f]['sep_oil_burn']))
    print()

    # ---- machine-readable comparison -------------------------------
    out = os.path.join(OUT_RESULTS, 'sensitivity_2024.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['feature', 'sep_oil_clean_all', 'sep_oil_clean_no2024',
                    'sep_oil_burn_all', 'sep_oil_burn_no2024', 'delta_burn',
                    'per_event_median_all', 'per_event_median_no2024',
                    'ci_lo_no2024', 'ci_hi_no2024'])
        for feat in full['rank_burn']:
            a = full['table'][feat]
            b = excl['table'].get(feat)
            if not b:
                continue
            ci = excl['ci'].get(feat, (np.nan, np.nan))
            w.writerow([feat,
                        '%.4f' % a['sep_oil_clean'], '%.4f' % b['sep_oil_clean'],
                        '%.4f' % a['sep_oil_burn'], '%.4f' % b['sep_oil_burn'],
                        '%.4f' % (b['sep_oil_burn'] - a['sep_oil_burn']),
                        '%.4f' % a['per_event_median'],
                        '%.4f' % b['per_event_median'],
                        '%.4f' % ci[0], '%.4f' % ci[1]])
    print('written: %s' % out)


if __name__ == '__main__':
    main()
