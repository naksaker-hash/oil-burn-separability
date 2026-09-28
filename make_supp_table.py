#!/usr/bin/env python3
"""
Supplementary table: all 50 screened features, from results/separability.csv.

Emits LaTeX to stdout. Longtable, so it breaks across pages.

Features that are monotone transforms of another in the set are marked
with a dagger: they rank pixels identically to their partner and so have
identical AUC by construction. They are kept in the table because the
literature reports them separately, but they are not independent
evidence and the counts in the text are given on the distinct set.

Usage:
    python3 make_supp_table.py > supp_table_features.tex
"""

import csv
import os
import sys

RESULTS = os.path.join('results', 'separability.csv')

# Redundant member -> the feature it duplicates. Established by Spearman
# correlation on all 34141 rows; see TABLE3_AND_RANKING.md.
DUPLICATE_OF = {
    'logBAI': 'BAI',
    'NDWI': 'NDBI',
    'FMR': 'NDBI',
    'd_NDWI': 'd_NDBI',
}

# LaTeX for feature names. Anything not here is escaped mechanically.
def texname(f):
    if f.startswith('d_'):
        return r'$\Delta$' + f[2:].replace('_', r'\_')
    return f.replace('_', r'\_')


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else RESULTS
    if not os.path.exists(path):
        raise SystemExit('Not found: %s. Run oil_burn_separability.py first.'
                         % path)
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    rows.sort(key=lambda r: -float(r['sep_oil_burn']))

    def num(r, k, nd=3):
        try:
            v = float(r[k])
        except (TypeError, ValueError):
            return '--'
        return '--' if v != v else '%.*f' % (nd, v)

    out = []
    out.append(r'\begin{longtable}{l r r r r r r r l}')
    out.append(r'\caption{All 50 screened features, ordered by '
               r'separability against burnt ground. Separability is '
               r'$|\mathrm{AUC}-0.5|\times 2$, where 0 means the two '
               r'classes are interchangeable on that feature. Per-event '
               r'columns '
               r'cover the 52 burn events contributing at least 100 '
               r'pixels. Intervals are 95\,\% percentile intervals from '
               r'a cluster bootstrap over burn events, computed for the '
               r'features ranking in the top twelve on either '
               r'comparison and for those named in the text. '
               r'JM is the Jeffries--Matusita distance on the same two '
               r'comparisons, tabulated for continuity with the remote '
               r'sensing literature. It assumes normally distributed '
               r'classes and no claim in the paper rests on it. '
               r'$\dagger$ marks a feature that is a '
               r'monotone transform of another in the set and therefore '
               r'has identical AUC by construction. The counts quoted '
               r'in the text are computed on the 46 distinct features.}'
               r'\label{tab:allfeatures}\\')
    out.append(r'\hline')
    hdr = (r'Feature & vs clean & vs burn & Drop & JM clean & JM burn & '
           r'Per-event & $>$0.80 & 95\,\% CI (burn) \\')
    out.append(hdr)
    out.append(r'\hline')
    out.append(r'\endfirsthead')
    out.append(r'\hline')
    out.append(hdr)
    out.append(r'\hline')
    out.append(r'\endhead')
    out.append(r'\hline')
    out.append(r'\endfoot')

    for r in rows:
        f = r['feature']
        name = texname(f)
        if f in DUPLICATE_OF:
            name += r'$^{\dagger}$'
        lo, hi = r['burn_ci_lo_event'], r['burn_ci_hi_event']
        try:
            ci = ('[%.3f, %.3f]' % (float(lo), float(hi))
                  if float(lo) == float(lo) else '')
        except (TypeError, ValueError):
            ci = ''
        out.append('%s & %s & %s & %s & %s & %s & %s & %s/%s & %s \\\\'
                   % (name, num(r, 'sep_oil_clean'), num(r, 'sep_oil_burn'),
                      num(r, 'drop'), num(r, 'jm_oil_clean', 2),
                      num(r, 'jm_oil_burn', 2), num(r, 'per_event_median'),
                      r['per_event_n_gt_080'], r['per_event_n'], ci))

    out.append(r'\end{longtable}')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
