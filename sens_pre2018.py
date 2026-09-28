#!/usr/bin/env python3
"""
Pre-2018 sensitivity run.

EFFIS perimeters are delineated from MODIS at 250 m before 2018 and
refined with Sentinel-2 at 20 m from 2018 onwards. The earlier events
are therefore less precisely bounded than the later ones, and the
limitations section says so without testing it. This drops every event
that burned before 2018 and re-runs the headline analysis.

The statistics come from oil_burn_separability.py and the comparison
machinery from sens_2024.py, both unchanged, so this run and the 2024
run differ only in which events they exclude.

Usage:
    python3 sens_pre2018.py [path/to/oil_burn_samples.csv]
"""

import csv
import os
import sys

import numpy as np

import oil_burn_separability as ob
from sens_2024 import analyse, year, fmt, HEADLINE, OUT_RESULTS

FIRST_REFINED_YEAR = 2018


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

    keep = lambda r: int(year(r)) >= FIRST_REFINED_YEAR
    burn_x = [r for r in burn if keep(r)]
    ctrl_x = [r for r in ctrl if keep(r)]

    print('=' * 70)
    print('PRE-%d SENSITIVITY RUN' % FIRST_REFINED_YEAR)
    print('=' * 70)
    print('dropped: %d burn px, %d control px, %d events (MODIS-era perimeters)'
          % (len(burn) - len(burn_x), len(ctrl) - len(ctrl_x),
             len({r['event_id'] for r in burn})
             - len({r['event_id'] for r in burn_x})))
    print()

    full = analyse(oil, edge, clean, burn, ctrl, feats, 'all events')
    print('run 1 of 2 done')
    excl = analyse(oil, edge, clean, burn_x, ctrl_x, feats, 'Sentinel-2 era only')
    print('run 2 of 2 done')
    print()

    print('%-30s %12s %14s' % ('', 'all events', '2018 onwards'))
    print('-' * 58)
    for k, lab in [('n_events', 'burn events'),
                   ('n_burn_px', 'burn pixels'),
                   ('n_big_events', 'events >= 100 px'),
                   ('n_gt090_clean', 'features > 0.90 vs clean'),
                   ('n_gt090_burn', 'features > 0.90 vs burn')]:
        print('%-30s %12d %14d' % (lab, full[k], excl[k]))
    for k, lab in [('median_clean', 'median sep. vs clean soil'),
                   ('median_burn', 'median sep. vs burn'),
                   ('qc_median', 'label quality, median AUC')]:
        print('%-30s %12s %14s' % (lab, fmt(full[k]), fmt(excl[k])))
    print()

    print('HEADLINE FEATURES, separability vs burn')
    print('%-10s %9s %12s %9s   %s' %
          ('feature', 'all', '2018 on', 'change', '95% CI, 2018 on'))
    print('-' * 62)
    for f in HEADLINE:
        a, b = full['table'].get(f), excl['table'].get(f)
        if not a or not b:
            continue
        ci = excl['ci'].get(f)
        cis = ('[%s, %s]' % (fmt(ci[0]), fmt(ci[1]))) if ci else ''
        print('%-10s %9.3f %12.3f %+9.3f   %s'
              % (f, a['sep_oil_burn'], b['sep_oil_burn'],
                 b['sep_oil_burn'] - a['sep_oil_burn'], cis))
    print()

    print('SUBSTRATE-MATCHED BURN PIXELS')
    print('  matched: %d (all), %d (2018 on)' % (full['n_matched'], excl['n_matched']))
    for f in ['MIRBI', 'NBR2', 'NDVI']:
        if f in full['matched'] and f in excl['matched']:
            print('  %-7s all %s -> %s   2018 on %s -> %s'
                  % (f, fmt(full['table'][f]['sep_oil_burn']), fmt(full['matched'][f]),
                     fmt(excl['table'][f]['sep_oil_burn']), fmt(excl['matched'][f])))
    print()

    print('TWO-FEATURE DISCRIMINANT, leave-one-burn-event-out')
    print('  pair            %s+%s' % (excl['det'], excl['dis']))
    print('  folds           %d (all) / %d (2018 on)' % (full['loeo_n'], excl['loeo_n']))
    print('  median          %s (all) / %s (2018 on)'
          % (fmt(full['loeo_median']), fmt(excl['loeo_median'])))
    print('  p10-p90         %s-%s / %s-%s'
          % (fmt(full['loeo_p10']), fmt(full['loeo_p90']),
             fmt(excl['loeo_p10']), fmt(excl['loeo_p90'])))
    print()

    out = os.path.join(OUT_RESULTS, 'sensitivity_pre2018.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['feature', 'sep_oil_burn_all', 'sep_oil_burn_2018on',
                    'delta', 'ci_lo_2018on', 'ci_hi_2018on'])
        for feat in full['rank_burn']:
            a = full['table'][feat]
            b = excl['table'].get(feat)
            if not b:
                continue
            ci = excl['ci'].get(feat, (np.nan, np.nan))
            w.writerow([feat, '%.4f' % a['sep_oil_burn'], '%.4f' % b['sep_oil_burn'],
                        '%.4f' % (b['sep_oil_burn'] - a['sep_oil_burn']),
                        '%.4f' % ci[0], '%.4f' % ci[1]])
    print('written: %s' % out)


if __name__ == '__main__':
    main()
