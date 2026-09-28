#!/usr/bin/env python3
"""
Label quality for a burn perimeter source, and the clean-soil ranking.

Two jobs, both feeding Table 3 and section 4.1.

1. LABEL QUALITY. For each event, how well its burn pixels separate from
   its OWN control ring on d_NBR. This says whether the perimeters
   actually contain burnt ground. It is a validation statistic, not a
   result: d_NBR is also a feature under test, which is legitimate here
   only because no event is ever excluded on it. Screening by it would
   be circular; measuring it is not.

   This is the statistic that showed MCD64A1 at 500 m to be inadequate
   (median 0.639, with events scoring BELOW 0.5 because the char sat in
   the control annulus rather than the labelled polygon) and EFFIS to be
   usable.

   Point it at a MODIS-perimeter export to fill in the unverified
   MCD64A1 counts in Table 3; the EFFIS numbers come from the archived
   oil_burn_samples.csv.

2. CLEAN-SOIL RANKING. The top features separating oil from clean soil,
   printed for comparison against paper 1's index screen (section 4.1
   claims this study's ranking reproduces it). The comparison itself has
   to be made by eye against paper 1's table.

Usage:
    python3 label_quality.py [path/to/samples.csv]
"""

import os
import sys
import numpy as np

import oil_burn_separability as ob


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else ob.CSV_DEFAULT
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)

    rows, feats = ob.load(path)
    burn = ob.subset(rows, 'burn')
    ctrl = ob.subset(rows, 'control_unburned')
    oil = ob.subset(rows, 'oil')
    clean = ob.subset(rows, 'control_clean')

    # Perimeter source is not recorded per row, so report what is here
    # and let the caller label it.
    ev = sorted({r['event_id'] for r in burn})
    print('=' * 60)
    print('LABEL QUALITY  burn vs its own control ring, per event, d_NBR')
    print('=' * 60)
    print('source file : %s' % path)
    print('events      : %d' % len(ev))

    scores = []
    for e in ev:
        b = [r for r in burn if r['event_id'] == e]
        c = [r for r in ctrl if r['event_id'] == e]
        a = ob.auc(ob.column(c, 'd_NBR'), ob.column(b, 'd_NBR'))
        if np.isfinite(a):
            scores.append(a)
    s = np.array(scores)
    if not s.size:
        raise SystemExit('No event scored. Are burn and control_unburned '
                         'both present, with matching event_id values?')

    print('events scored : %d' % s.size)
    print('median AUC    : %.3f' % np.median(s))
    print('mean AUC      : %.3f' % s.mean())
    print('  >= 0.90     : %d of %d' % ((s >= 0.90).sum(), s.size))
    print('  >  0.80     : %d of %d' % ((s > 0.80).sum(), s.size))
    print('  <  0.50     : %d of %d' % ((s < 0.50).sum(), s.size))
    print('min, max      : %.3f, %.3f' % (s.min(), s.max()))
    print()

    print('=' * 60)
    print('TOP 15, oil vs CLEAN SOIL  (for comparison with paper 1)')
    print('=' * 60)
    oil_c = {f: ob.column(oil, f) for f in feats}
    clean_c = {f: ob.column(clean, f) for f in feats}
    rank = []
    for f in feats:
        v = ob.separability(oil_c[f], clean_c[f])
        if np.isfinite(v):
            rank.append((v, f))
    rank.sort(reverse=True)
    print('%-4s %-12s %10s' % ('rank', 'feature', 'separability'))
    print('-' * 30)
    for i, (v, f) in enumerate(rank[:15], start=1):
        print('%-4d %-12s %10.4f' % (i, f, v))
    print()
    print('%d of %d features exceed 0.90 against clean soil.'
          % (sum(1 for v, _ in rank if v > 0.90), len(rank)))


if __name__ == '__main__':
    main()
