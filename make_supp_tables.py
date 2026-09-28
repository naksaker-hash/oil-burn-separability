#!/usr/bin/env python3
"""
Supplementary tables S3-S6, emitted as one LaTeX fragment.

  S3  per-event acquisition intervals and separability
  S4  oil-mask sensitivity
  S5  nested leave-one-fire-out validation and its baselines
  S6  the full matched set against the 2018-onward subset

Inputs are the CSVs written by temporal_lags.py, oil_mask_sensitivity.py,
nested_loeo.py and sens_pre2018.py, plus the sample table itself for the
event-level descriptives in S6.

Usage:
    python3 make_supp_tables.py > overleaf/supp_tables_extra.tex
"""

import csv
import os
import numpy as np

from oil_burn_separability import (load, subset, column, separability,
                                   auc, MIN_EVENT_N)

R = 'results'
OUT = []


def w(s=''):
    OUT.append(s)


def f(x, nd=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return '--'
    return '--' if v != v else '%.*f' % (nd, v)


def texname(name):
    if name.startswith('d_'):
        return r'$\Delta$' + name[2:].replace('logBAI', r'$\log$\,BAI') \
            if name[2:] != 'logBAI' else r'$\Delta\log$\,BAI'
    return name


# ---------------------------------------------------------------- S3
def table_s3():
    rows = list(csv.DictReader(open(os.path.join(R, 'event_intervals.csv'))))
    rows.sort(key=lambda r: -float(r['dt_tot']))
    w(r'\begin{longtable}{l r r r r r r r}')
    w(r'\caption{Acquisition intervals and separability from oil, for '
      r'each of the 52 burn events contributing at least 100 sampled '
      r'pixels. $\Delta t_{\mathrm{pre}}$ is the number of days from '
      r'the pre-event scene to the start of the fire, '
      r'$\Delta t_{\mathrm{post}}$ the number from the end of the fire '
      r'to the post-event scene, and $\Delta t_{\mathrm{tot}}$ the '
      r'separation between the two scenes. Rows are ordered by '
      r'$\Delta t_{\mathrm{tot}}$. The oil site, for comparison, has '
      r'$\Delta t_{\mathrm{pre}} = 4$, $\Delta t_{\mathrm{post}} = 1$ '
      r'and $\Delta t_{\mathrm{tot}} = 5$ days. Rank correlations '
      r'between $\Delta t_{\mathrm{tot}}$ and the separability columns '
      r'are $+0.17$ ($\Delta$NBR), $-0.06$ ($\Delta\log$\,BAI), '
      r'$-0.04$ (MIRBI) and $-0.05$ (NBR2).}'
      r'\label{tab:intervals}\\')
    hdr = (r'Event & Pixels & Area (ha) & $\Delta t_{\mathrm{pre}}$ & '
           r'$\Delta t_{\mathrm{post}}$ & $\Delta t_{\mathrm{tot}}$ & '
           r'$S$ MIRBI & $S$ $\Delta$NBR \\')
    for tag in (r'\endfirsthead', r'\endhead'):
        w(r'\hline'); w(hdr); w(r'\hline'); w(tag)
    w(r'\hline'); w(r'\endfoot')
    for r in rows:
        w(r'\texttt{%s} & %s & %s & %d & %d & %d & %s & %s \\' % (
            r['event_id'].replace('_', r'\_'), r['n_pixels'],
            f(r['area_ha'], 0), int(float(r['dt_pre'])),
            int(float(r['dt_post'])), int(float(r['dt_tot'])),
            f(r['S_MIRBI']), f(r['S_d_NBR'])))
    w(r'\end{longtable}')
    w()


# ---------------------------------------------------------------- S4
def table_s4():
    rows = list(csv.DictReader(open(os.path.join(R,
                                                 'oil_mask_sensitivity.csv'))))
    w(r'\begin{table}[htbp]')
    w(r'  \centering')
    w(r'  \caption{Sensitivity of the main comparisons to how the plume '
      r'boundary was drawn. The core mask keeps the 164 sampled oil '
      r'pixels whose four neighbours are all themselves within the '
      r'plume (6.56~ha). The plume mask is the 218 pixels used '
      r'throughout the paper (8.72~ha). The inclusive mask adds the 100 '
      r'plume-edge pixels (12.72~ha). The masks are nested. '
      r'Discrimination against burn scars varies by at most 0.04 across '
      r'the three, well inside the cluster-bootstrap interval, while '
      r'detection against clean soil weakens as partly contaminated '
      r'margin pixels are admitted, as expected.}')
    w(r'  \label{tab:masks}')
    w(r'  \begin{tabular}{lcccccc}')
    w(r'    \toprule')
    w(r'    & \multicolumn{3}{c}{vs.\ clean soil} & '
      r'\multicolumn{3}{c}{vs.\ burn scar} \\')
    w(r'    \cmidrule(lr){2-4}\cmidrule(lr){5-7}')
    w(r'    Feature & Core & Plume & Inclusive & Core & Plume & Inclusive \\')
    w(r'    \midrule')
    for r in rows:
        w(r'    %s & %s & %s & %s & %s & %s & %s \\' % (
            texname(r['feature']), f(r['core_clean']), f(r['plume_clean']),
            f(r['inclusive_clean']), f(r['core_burn']), f(r['plume_burn']),
            f(r['inclusive_burn'])))
    w(r'    \bottomrule')
    w(r'  \end{tabular}')
    w(r'\end{table}')
    w()


# ---------------------------------------------------------------- S5
def table_s5():
    rows = list(csv.DictReader(open(os.path.join(R, 'nested_loeo.csv'))))
    w(r'\begin{table}[htbp]')
    w(r'  \centering')
    w(r'  \caption{Nested leave-one-fire-out validation. For each '
      r'withheld fire the detector and the discriminator are selected '
      r'using only the remaining fires and the fitting half of the oil '
      r'pixels, the discriminant is fitted on those two features, and '
      r'the withheld fire is scored against the withheld oil pixels. '
      r'The 52 fire-level test units are each evaluated under five '
      r'random splits of the oil pixels, giving 260 evaluations. The '
      r'repetitions measure sensitivity to the oil split and do not add '
      r'independent validation. The fixed pair rows repeat the '
      r'procedure with $\Delta$NBR and MIRBI imposed rather than '
      r'selected. They are identical because the nested selection chose '
      r'that pair in all 260 evaluations.}')
    w(r'  \label{tab:nested}')
    w(r'  \begin{tabular}{lccc}')
    w(r'    \toprule')
    w(r'    Model & Median $S$ & 10th pct. & 90th pct. \\')
    w(r'    \midrule')
    spec = [('nested_pair', 'Nested pair, selected in fold'),
            ('nested_dis', r'\quad its discriminator alone'),
            ('nested_det', r'\quad its detector alone'),
            ('fixed_pair', r'Fixed pair, $\Delta$NBR $+$ MIRBI'),
            ('fixed_dis', r'\quad MIRBI alone'),
            ('fixed_det', r'\quad $\Delta$NBR alone')]
    for k, lab in spec:
        a = np.array([float(r[k]) for r in rows if r[k] not in ('', None)])
        a = a[np.isfinite(a)]
        w(r'    %s & %.3f & %.3f & %.3f \\' % (
            lab, np.median(a), np.percentile(a, 10), np.percentile(a, 90)))
    w(r'    \bottomrule')
    w(r'  \end{tabular}')
    w(r'\end{table}')
    w()


# ---------------------------------------------------------------- S6
def table_s6():
    rows, feats = load('oil_burn_samples.csv')
    oil = subset(rows, 'oil')
    burn = subset(rows, 'burn')
    ctrl = subset(rows, 'control_unburned')

    def year(r):
        return int(str(r['event_date'])[:4])

    def stats(keep):
        b = [r for r in burn if keep(r)]
        ev = sorted({r['event_id'] for r in b})
        big = {e: [r for r in b if r['event_id'] == e] for e in ev}
        big = {e: v for e, v in big.items() if len(v) >= MIN_EVENT_N}
        area = np.array([float(big[e][0]['area_ha']) for e in big])
        ndvi = np.array([float(big[e][0]['pre_ndvi']) for e in big])
        slope = np.array([float(big[e][0]['slope_mean']) for e in big])
        qc = []
        for e in ev:
            bb = [r for r in b if r['event_id'] == e]
            cc = [r for r in ctrl if r['event_id'] == e and keep(r)]
            if len(bb) >= 10 and len(cc) >= 10:
                a = auc(column(cc, 'd_NBR'), column(bb, 'd_NBR'))
                if np.isfinite(a):
                    qc.append(a)
        out = {
            'events': len(ev), 'big': len(big), 'pixels': len(b),
            'area': np.median(area), 'ndvi': np.median(ndvi),
            'slope': np.median(slope), 'qc': np.median(qc) if qc else np.nan,
        }
        for f_ in ('MIRBI', 'NBR2', 'd_NBR', 'd_logBAI'):
            out[f_] = separability(column(oil, f_), column(b, f_))
            per = [separability(column(v, f_), column(oil, f_))
                   for v in big.values()]
            per = [x for x in per if np.isfinite(x)]
            out['pe_' + f_] = np.median(per) if per else np.nan
        return out

    a = stats(lambda r: True)
    c = stats(lambda r: year(r) >= 2018)

    w(r'\begin{table}[htbp]')
    w(r'  \centering')
    w(r'  \caption{The full matched set against the subset from 2018 '
      r'onwards, for which EFFIS perimeters are refined with '
      r'Sentinel-2 at 20~m rather than mapped from MODIS at 250~m. '
      r'Label quality is the median per-event AUC of burn pixels '
      r'against its own control annulus. '
      r'Per-event medians cover the events contributing at least 100 '
      r'pixels. The subset is better labelled and separates further '
      r'from oil on every discriminating feature, so the full-sample '
      r'values reported in the paper are the conservative ones.}')
    w(r'  \label{tab:pre2018}')
    w(r'  \begin{tabular}{lcc}')
    w(r'    \toprule')
    w(r'    & All events & 2018 onwards \\')
    w(r'    \midrule')
    rowspec = [
        ('Events', '%d', 'events'), ('Events with $\\geq$100 pixels', '%d', 'big'),
        ('Burn pixels', '%d', 'pixels'),
        ('Median area (ha)', '%.0f', 'area'),
        ('Median pre-fire NDVI', '%.3f', 'ndvi'),
        ('Median slope (deg)', '%.2f', 'slope'),
        ('Label quality (median AUC)', '%.3f', 'qc'),
        (r'\midrule', None, None),
        ('MIRBI, pooled', '%.3f', 'MIRBI'),
        ('MIRBI, per-event median', '%.3f', 'pe_MIRBI'),
        ('NBR2, pooled', '%.3f', 'NBR2'),
        ('NBR2, per-event median', '%.3f', 'pe_NBR2'),
        (r'$\Delta$NBR, pooled', '%.3f', 'd_NBR'),
        (r'$\Delta\log$\,BAI, pooled', '%.3f', 'd_logBAI'),
    ]
    for lab, fmt, key in rowspec:
        if fmt is None:
            w('    ' + lab)
            continue
        w(('    %s & ' + fmt + ' & ' + fmt + r' \\') % (lab, a[key], c[key]))
    w(r'    \bottomrule')
    w(r'  \end{tabular}')
    w(r'\end{table}')


# ---------------------------------------------------------------- S7
def table_s7():
    rows = list(csv.DictReader(open(os.path.join(R,
                                                 'sequential_screen.csv'))))
    P = 0.0097

    def med(k):
        a = np.array([float(r[k]) for r in rows if r[k] not in ('', None)])
        a = a[np.isfinite(a)]
        return np.median(a) if a.size else float('nan')

    w(r'\begin{table}[htbp]')
    w(r'  \centering')
    w(r'  \caption{The two-stage screen of Sect.~3.7 of the main text. '
      r'Stage~1 thresholds $\Delta$NBR to retain 90\,\% of the training '
      r'oil pixels, and stage~2 thresholds MIRBI to retain 90\,\% of the '
      r'stage~1 survivors. Both thresholds are fitted with the test '
      r'fire, its own control annulus and half the oil pixels withheld, '
      r'and are quantiles of the training oil half rather than values '
      r'tuned against the outcome. The upper block gives the median '
      r'fraction of each withheld class that survives. The lower block '
      r'substitutes those rates into Eq.~(1) of the main text at an '
      r'assumed scene composition. Because that composition is assumed '
      r'and no scene is classified, the precision figures illustrate the '
      r'arithmetic rather than measure performance.}')
    w(r'  \label{tab:sequential}')
    w(r'  \begin{tabular}{lcc}')
    w(r'    \toprule')
    w(r'    Withheld class & After stage 1 & After stage 2 \\')
    w(r'    \midrule')
    for lab, a, b in (('Oil', 'oil1', 'oil2'),
                      ('Burn, the withheld fire', 'burn1', 'burn2'),
                      ('Clean soil', 'clean1', 'clean2'),
                      ('Unburned control', 'ctrl1', 'ctrl2')):
        w(r'    %s & %.3f & %.3f \\' % (lab, med(a), med(b)))
    w(r'    \midrule')
    w(r'    \multicolumn{3}{l}{\emph{Illustrative precision at '
      r'$P = 0.0097$}} \\')
    w(r'    Burnt fraction of scene & One stage & Two stages \\')
    bg1 = 0.5 * (med('clean1') + med('ctrl1'))
    bg2 = 0.5 * (med('clean2') + med('ctrl2'))
    for bf in (0.01, 0.05, 0.10):
        rest = 1.0 - P - bf
        fp1 = med('burn1') * bf + bg1 * rest
        fp2 = med('burn2') * bf + bg2 * rest
        p1 = med('oil1') * P / (med('oil1') * P + fp1)
        p2 = med('oil2') * P / (med('oil2') * P + fp2)
        w(r'    %.0f\,\%% & %.3f & %.3f \\' % (100 * bf, p1, p2))
    w(r'    \bottomrule')
    w(r'  \end{tabular}')
    w(r'\end{table}')


def main():
    table_s3()
    w(r'\clearpage')
    table_s4()
    table_s5()
    w(r'\clearpage')
    table_s6()
    table_s7()
    print('\n'.join(OUT))


if __name__ == '__main__':
    main()
