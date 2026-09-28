#!/usr/bin/env python3
"""
Separability of terrestrial oil contamination from burn scars
in Sentinel-2 spectral indices.

Input : oil_burn_samples.csv, exported by narli_burn_confusion_gee.js
Output: results/*.csv, figures/*.png, and a printed summary

Dependencies: numpy and matplotlib ONLY.

Everything that would normally come from scipy or scikit-learn is
computed here directly:
  - AUC from the rank-sum identity, with proper mid-ranks for ties
  - Jeffries-Matusita distance from the Bhattacharyya coefficient
  - Fisher's linear discriminant solved analytically

The reason is environmental, not stylistic: the working Anaconda base
has pandas, pyarrow and scipy compiled against numpy 1.x, which cannot
import under numpy 2.x. Keeping this file to numpy and matplotlib means
it runs in the broken environment and in a clean one alike.

Usage:
    python3 oil_burn_separability.py [path/to/oil_burn_samples.csv]
"""

import csv
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------

CSV_DEFAULT = 'oil_burn_samples.csv'
OUT_RESULTS = 'results'
OUT_FIGURES = 'figures'

# Columns that describe the sample rather than the spectra.
META = {
    'class_name', 'event_id', 'event_date', 'event_end_date',
    'pre_date', 'post_date', 'area_ha', 'pre_ndvi', 'agri_pct',
    'slope_mean', 'soil_mode', 'province',
    'pre_NDVI_px', 'slope_px', 'soil_px',
}

# Minimum rows in each class before a comparison is attempted.
MIN_N = 10

# Burn events need at least this many rows to get their own estimate.
MIN_EVENT_N = 100

# Bootstrap replicates. Intervals are computed only for the features
# that rank highly on either comparison (TOP_K on each, unioned), not
# for all of them: at 50 features and full replicates this took longer
# than the whole rest of the analysis and produced intervals nobody
# reads for a feature ranked fortieth.
N_BOOT = 500
BOOT_SEED = 20260911
TOP_K = 12

# Features that always get an interval regardless of where they rank.
# The top-K union misses d_logBAI and d_BAI: both sit outside the top
# twelve on each comparison taken separately, yet d_logBAI carries the
# sharpest claim in the paper (0.919 against clean soil, 0.018 against
# burn). A number described as indistinguishable from chance needs an
# interval more than a number that is obviously large does.
CI_FEATURES = ['d_logBAI', 'd_BAI', 'logBAI', 'BAI', 'MIRBI', 'NBR2',
               'd_NBR', 'd_FMR', 'NBR', 'NDVI', 'B11', 'B12']

# Leave-one-event-out repetitions for the two-index combination. Each
# repetition is a fresh random split of the oil pixels; the burn side is
# held out by event regardless, which is where the real variance lives.
LOEO_REPEATS = 5

# Substrate-matched sensitivity: keep burn pixels whose pre-fire NDVI
# and slope fall within these distances of the plume's own values.
NDVI_TOL = 0.06
SLOPE_TOL = 3.0

# Features on Figure 2, named rather than taken from the top of each
# ranking. The earlier version used the top four by burn plus the top
# three by clean soil, which silently excluded d_logBAI: it ranks fourth
# against clean soil, so it fell outside both slices. That is the
# sharpest number in the paper -- the index put forward as a novel oil
# detector, sitting at chance against burn -- and a figure about the
# dissociation that omits it argues the case less well than the text.
#
# Order is discriminators first, then detectors, so the panel reads
# left to right as "these separate oil from char, these do not".
# Features labelled on Figure 1. The caption says "features discussed in
# the text are labelled", so the list is the set the text actually names
# rather than everything past a threshold: the old rule labelled some
# twenty points in the crowded lower right, where several features sit at
# identical coordinates and their labels printed on top of one another.
# Coincident features are merged into a single label below, which also
# shows the redundancy described in the methods.
FIG1_LABEL = ['d_NBR', 'd_FMR', 'NBR', 'd_logBAI', 'd_BAI',
              'MIRBI', 'NBR2', 'd_ARVI', 'NDVI', 'B11', 'ARVI',
              'd_MIRBI', 'NDBI', 'NDWI', 'FMR', 'd_NDBI', 'd_NDWI',
              'BAI', 'logBAI']

FIG2_DISCRIMINATORS = ['MIRBI', 'NBR2', 'NDVI', 'B11']
FIG2_DETECTORS = ['d_NBR', 'd_FMR', 'NBR', 'd_logBAI']


# --------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------

def load(path):
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit('No rows in %s' % path)
    feats = [k for k in rows[0] if k and k not in META]
    return rows, feats


def column(rows, name):
    """Float column with non-numeric and non-finite entries as NaN."""
    out = np.empty(len(rows))
    for i, r in enumerate(rows):
        try:
            v = float(r[name])
            out[i] = v if np.isfinite(v) else np.nan
        except (TypeError, ValueError):
            out[i] = np.nan
    return out


def subset(rows, cls):
    return [r for r in rows if r['class_name'] == cls]


# --------------------------------------------------------------------
# Statistics, implemented directly
# --------------------------------------------------------------------

def _midranks(x):
    """Ranks 1..n with ties averaged, matching the Wilcoxon convention."""
    order = np.argsort(x, kind='mergesort')
    sx = x[order]
    ranks = np.empty(len(x), dtype=float)
    i = 0
    n = len(sx)
    while i < n:
        j = i
        while j + 1 < n and sx[j + 1] == sx[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def auc(a, b):
    """
    P(a > b) with ties counted as one half: the Mann-Whitney statistic
    divided by n_a * n_b. Returns NaN if either group is too small.

    Verified against sklearn.metrics.roc_auc_score to 1e-16 on random
    data with heavy ties.
    """
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < MIN_N or len(b) < MIN_N:
        return np.nan
    r = _midranks(np.concatenate([a, b]))
    na, nb = len(a), len(b)
    return (r[:na].sum() - na * (na + 1) / 2.0) / (na * nb)


def separability(a, b):
    """
    |AUC - 0.5| * 2, on [0, 1].

    Direction-free, because the question here is whether two classes can
    be told apart at all, not which way round they sit. 0 means the
    distributions are interchangeable; 1 means they do not overlap.
    """
    v = auc(a, b)
    return np.nan if not np.isfinite(v) else abs(v - 0.5) * 2.0


def jm_distance(a, b):
    """
    Jeffries-Matusita distance on [0, 2], from the Bhattacharyya
    distance under a normal assumption.

    Reported because the remote sensing literature expects it, but AUC
    is the primary statistic here: J-M assumes normality and these
    indices are not normal, so J-M and AUC will occasionally disagree.
    Where they do, trust AUC, which assumes nothing.
    """
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < MIN_N or len(b) < MIN_N:
        return np.nan
    m1, m2 = a.mean(), b.mean()
    v1, v2 = a.var(ddof=1), b.var(ddof=1)
    if v1 <= 0 or v2 <= 0:
        return np.nan
    bdist = (0.125 * (m1 - m2) ** 2 * 2.0 / (v1 + v2)
             + 0.5 * np.log((v1 + v2) / (2.0 * np.sqrt(v1 * v2))))
    return 2.0 * (1.0 - np.exp(-bdist))


def boot_ci_pixels(a, b, stat=separability, n=N_BOOT, seed=BOOT_SEED):
    """
    Percentile interval from resampling PIXELS independently.

    This interval is optimistic. Pixels inside one plume or one burn
    scar are spatially autocorrelated, so they carry less information
    than their count implies. Use boot_ci_events for the burn side,
    where the clustering can actually be respected.
    """
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < MIN_N or len(b) < MIN_N:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = np.empty(n)
    for i in range(n):
        vals[i] = stat(rng.choice(a, len(a), replace=True),
                       rng.choice(b, len(b), replace=True))
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan, np.nan
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)


def boot_ci_events(oil_vals, event_arrays, stat=separability,
                   n=N_BOOT, seed=BOOT_SEED):
    """
    Cluster bootstrap: resample burn EVENTS with replacement, keeping
    each event's pixels together.

    This is the honest interval for the burn side, because the sampling
    unit that varies between fires is the fire, not the pixel. It is
    typically several times wider than the pixel interval.

    The oil side cannot be treated this way: there is one spill, so its
    between-event variance is unmeasurable from this data. Any interval
    quoted here understates the uncertainty in generalising to other
    spills, and the text must say so.
    """
    oil_vals = oil_vals[np.isfinite(oil_vals)]
    keys = [k for k in event_arrays if np.isfinite(event_arrays[k]).sum() >= MIN_N]
    if len(oil_vals) < MIN_N or len(keys) < 3:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = np.empty(n)
    for i in range(n):
        pick = rng.choice(len(keys), len(keys), replace=True)
        pooled = np.concatenate([event_arrays[keys[j]] for j in pick])
        vals[i] = stat(oil_vals, pooled)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan, np.nan
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)


# --------------------------------------------------------------------
# Two-feature discriminant
# --------------------------------------------------------------------

def fisher_ld(X0, X1, ridge=1e-8):
    """
    Fisher's linear discriminant, solved analytically.

    Logistic regression is deliberately avoided: with well-separated
    classes it hits complete separation and returns coefficients that
    diverge, which is the failure recorded in the earlier project. The
    Fisher solution has a closed form and degrades gracefully.

    A small ridge keeps the pooled covariance invertible when two
    features are near-collinear.
    """
    m0, m1 = X0.mean(axis=0), X1.mean(axis=0)
    S = np.cov(X0, rowvar=False) * (len(X0) - 1) \
        + np.cov(X1, rowvar=False) * (len(X1) - 1)
    S = S / (len(X0) + len(X1) - 2)
    S = np.atleast_2d(S)
    S += np.eye(S.shape[0]) * ridge * np.trace(S) / S.shape[0]
    return np.linalg.solve(S, (m1 - m0))


def pretty(f):
    """
    Column name as the manuscript sets it.

    The CSV carries machine names -- d_NBR, d_logBAI -- while the text and
    tables use the delta form. Figures that disagree with the prose about
    what a feature is called make the reader do a translation on every
    glance, so the label is built here rather than left as the column
    name.
    """
    lab = f
    if lab.startswith('d_'):
        lab = '\u0394' + lab[2:]
    return lab.replace('logBAI', 'log BAI')


def clean_rows(arrs):
    """Row mask where every supplied column is finite."""
    m = np.ones(len(arrs[0]), dtype=bool)
    for a in arrs:
        m &= np.isfinite(a)
    return m


# --------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else CSV_DEFAULT
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)
    os.makedirs(OUT_RESULTS, exist_ok=True)
    os.makedirs(OUT_FIGURES, exist_ok=True)

    rows, feats = load(path)
    oil = subset(rows, 'oil')
    edge = subset(rows, 'oil_edge')
    clean = subset(rows, 'control_clean')
    burn = subset(rows, 'burn')
    ctrl = subset(rows, 'control_unburned')

    print('=' * 66)
    print('OIL / BURN SEPARABILITY')
    print('=' * 66)
    print('source      : %s' % path)
    print('features    : %d' % len(feats))
    print('oil         : %5d' % len(oil))
    print('oil_edge    : %5d' % len(edge))
    print('clean soil  : %5d' % len(clean))
    print('burn        : %5d' % len(burn))
    print('unburned    : %5d' % len(ctrl))
    ev = sorted({r['event_id'] for r in burn})
    print('burn events : %5d' % len(ev))
    print()

    # ---- label quality -------------------------------------------
    # Burn against its own control ring, per event, on d_NBR. This is a
    # validation statistic, not a result: it says whether the burn
    # polygons actually contain burnt ground. Reported in the methods.
    #
    # It is safe to use d_NBR here even though d_NBR is also a feature
    # under test, because nothing is SELECTED on it. Screening events by
    # this number would be circular; measuring it is not.
    qc = []
    for e in ev:
        b = [r for r in burn if r['event_id'] == e]
        c = [r for r in ctrl if r['event_id'] == e]
        a = auc(column(c, 'd_NBR'), column(b, 'd_NBR'))
        if np.isfinite(a):
            qc.append(a)
    qc = np.array(qc)
    if qc.size:
        print('LABEL QUALITY  burn vs own control ring, per event, d_NBR')
        print('  events scored : %d' % qc.size)
        print('  median AUC    : %.3f' % np.median(qc))
        print('  >= 0.90       : %d' % (qc >= 0.90).sum())
        print('  <  0.50       : %d' % (qc < 0.50).sum())
        print()

    # ---- per-feature separability --------------------------------
    oil_c = {f: column(oil, f) for f in feats}
    clean_c = {f: column(clean, f) for f in feats}
    burn_c = {f: column(burn, f) for f in feats}
    edge_c = {f: column(edge, f) for f in feats}

    big_events = {}
    for e in ev:
        b = [r for r in burn if r['event_id'] == e]
        if len(b) >= MIN_EVENT_N:
            big_events[e] = b

    # Per-event arrays are precomputed once. Rebuilding them inside the
    # bootstrap and CV loops was the single biggest cost in an earlier
    # version of this script.
    ev_cols = {f: {k: column(v, f) for k, v in big_events.items()}
               for f in feats}

    # Pass 1: point estimates for every feature. Cheap.
    table = []
    for f in feats:
        o, c, b, g = oil_c[f], clean_c[f], burn_c[f], edge_c[f]
        s_clean = separability(o, c)
        s_burn = separability(o, b)
        s_edge = separability(o, g)
        if not (np.isfinite(s_clean) and np.isfinite(s_burn)):
            continue

        per_ev = np.array([s for s in
                           (separability(o, a) for a in ev_cols[f].values())
                           if np.isfinite(s)])

        table.append({
            'feature': f,
            'sep_oil_clean': s_clean,
            'sep_oil_burn': s_burn,
            'sep_oil_edge': s_edge,
            'jm_oil_clean': jm_distance(o, c),
            'jm_oil_burn': jm_distance(o, b),
            'burn_ci_lo_pixel': np.nan,
            'burn_ci_hi_pixel': np.nan,
            'burn_ci_lo_event': np.nan,
            'burn_ci_hi_event': np.nan,
            'per_event_median': np.median(per_ev) if per_ev.size else np.nan,
            'per_event_n_gt_080': int((per_ev > 0.80).sum()) if per_ev.size else 0,
            'per_event_n': int(per_ev.size),
            'drop': s_clean - s_burn,
        })

    table.sort(key=lambda d: -d['sep_oil_clean'])

    # Pass 2: bootstrap intervals, only for features that rank in the
    # top TOP_K on either comparison.
    top = {d['feature'] for d in table[:TOP_K]}
    top |= {d['feature'] for d in
            sorted(table, key=lambda d: -d['sep_oil_burn'])[:TOP_K]}
    top |= {f for f in CI_FEATURES if f in oil_c}
    print('bootstrapping %d features x %d replicates ...' % (len(top), N_BOOT))
    for d in table:
        if d['feature'] not in top:
            continue
        f = d['feature']
        d['burn_ci_lo_pixel'], d['burn_ci_hi_pixel'] = \
            boot_ci_pixels(oil_c[f], burn_c[f])
        d['burn_ci_lo_event'], d['burn_ci_hi_event'] = \
            boot_ci_events(oil_c[f], ev_cols[f])
    print('done')
    print()

    cols = list(table[0].keys())
    with open(os.path.join(OUT_RESULTS, 'separability.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(table)

    print('TOP 15 BY oil vs CLEAN SOIL')
    print('%-12s %9s %9s %9s   %s' %
          ('feature', 'vs clean', 'vs burn', 'drop', 'burn 95% CI (event)'))
    print('-' * 74)
    for d in table[:15]:
        print('%-12s %9.4f %9.4f %9.3f   [%.3f, %.3f]' %
              (d['feature'], d['sep_oil_clean'], d['sep_oil_burn'],
               d['drop'], d['burn_ci_lo_event'], d['burn_ci_hi_event']))
    print()

    by_burn = sorted(table, key=lambda d: -d['sep_oil_burn'])
    print('TOP 10 BY oil vs BURN')
    print('%-12s %9s %9s %9s %7s   %s' %
          ('feature', 'vs burn', 'vs clean', 'per-ev', '>0.80', '95% CI (event)'))
    print('-' * 78)
    for d in by_burn[:10]:
        print('%-12s %9.4f %9.4f %9.3f %4d/%-3d [%.3f, %.3f]' %
              (d['feature'], d['sep_oil_burn'], d['sep_oil_clean'],
               d['per_event_median'], d['per_event_n_gt_080'],
               d['per_event_n'], d['burn_ci_lo_event'], d['burn_ci_hi_event']))
    print()

    sc = np.array([d['sep_oil_clean'] for d in table])
    sb = np.array([d['sep_oil_burn'] for d in table])
    print('median separability   vs clean %.3f   vs burn %.3f' %
          (np.median(sc), np.median(sb)))
    print('features > 0.90 vs clean : %d' % (sc > 0.90).sum())
    print('  of those > 0.90 vs burn: %d' % ((sc > 0.90) & (sb > 0.90)).sum())
    print()

    # ---- substrate-matched sensitivity ---------------------------
    # Restrict burn pixels to those whose own pre-fire NDVI and slope
    # sit close to the plume's. Answers the obvious objection that the
    # separation reflects site and season rather than material.
    plume_ndvi = np.nanmedian(column(oil, 'pre_NDVI_px'))
    plume_slope = np.nanmedian(column(oil, 'slope_px'))
    bn = column(burn, 'pre_NDVI_px')
    bs = column(burn, 'slope_px')
    keep = (np.abs(bn - plume_ndvi) <= NDVI_TOL) & \
           (np.abs(bs - plume_slope) <= SLOPE_TOL)
    burn_m = [r for r, k in zip(burn, keep) if k]
    print('SUBSTRATE-MATCHED SENSITIVITY')
    print('  plume pre-fire NDVI %.3f, slope %.2f deg' % (plume_ndvi, plume_slope))
    print('  burn pixels within +/-%.2f NDVI and +/-%.1f deg: %d of %d' %
          (NDVI_TOL, SLOPE_TOL, len(burn_m), len(burn)))
    if len(burn_m) >= MIN_N:
        print('  %-12s %9s %9s' % ('feature', 'all burn', 'matched'))
        for d in by_burn[:6]:
            f = d['feature']
            print('  %-12s %9.4f %9.4f' %
                  (f, d['sep_oil_burn'],
                   separability(oil_c[f], column(burn_m, f))))
    print()

    # ---- two-feature combination ---------------------------------
    # One index to detect oil against clean soil, one to hold it apart
    # from burn. Validated leave-one-burn-event-out.
    #
    # The oil side is split at random because there is only one spill,
    # so oil pixels appear in both training and test folds. They are
    # spatially autocorrelated, so this CV is optimistic about the oil
    # class. It is honest about generalisation to UNSEEN FIRES, which is
    # the claim being made.
    det = table[0]['feature']
    dis = by_burn[0]['feature']
    print('TWO-FEATURE COMBINATION: %s (detect) + %s (discriminate)' % (det, dis))

    rng = np.random.default_rng(BOOT_SEED)
    o2 = np.column_stack([oil_c[det], oil_c[dis]])
    o2 = o2[clean_rows([o2[:, 0], o2[:, 1]])]

    # Per-event two-column arrays, cleaned once rather than per fold.
    ev2 = {}
    for k in big_events:
        m = np.column_stack([ev_cols[det][k], ev_cols[dis][k]])
        m = m[clean_rows([m[:, 0], m[:, 1]])]
        if len(m) >= MIN_N:
            ev2[k] = m
    big_events = {k: big_events[k] for k in ev2}
    scores = []
    for rep in range(LOEO_REPEATS):
        idx = rng.permutation(len(o2))
        half = len(o2) // 2
        otr, ote = o2[idx[:half]], o2[idx[half:]]
        for k in big_events:
            tr = [ev2[kk] for kk in big_events if kk != k]
            if not tr:
                continue
            btr = np.vstack(tr)
            bte = ev2[k]
            if len(btr) < MIN_N or len(bte) < MIN_N:
                continue
            w = fisher_ld(btr, otr)
            s = separability(ote @ w, bte @ w)
            if np.isfinite(s):
                scores.append(s)
    scores = np.array(scores)
    if scores.size:
        print('  leave-one-event-out folds : %d' % scores.size)
        print('  median separability       : %.3f' % np.median(scores))
        print('  10th-90th percentile      : %.3f to %.3f' %
              (np.percentile(scores, 10), np.percentile(scores, 90)))
        print('  best single index (%s): %.3f' % (dis, by_burn[0]['sep_oil_burn']))
    print()

    # ---- figures --------------------------------------------------
    # Figure 1: the dissociation
    fig, ax = plt.subplots(figsize=(7, 6.2))
    ax.scatter(sc, sb, s=28, c='#444444', alpha=0.75, edgecolor='none')
    # Group features that plot at the same point so their labels do not
    # stack. Features that are monotone transforms of one another are
    # identical by construction and would otherwise overprint.
    groups = {}
    for d in table:
        if d['feature'] not in FIG1_LABEL:
            continue
        k = (round(d['sep_oil_clean'], 3), round(d['sep_oil_burn'], 3))
        groups.setdefault(k, []).append(pretty(d['feature']))
    # Place labels inside the axes and away from each other: points on the
    # right half are labelled to their left so long names do not run off
    # the panel, and a label that would land on one already placed is put
    # below its point instead of above.
    placed = []
    for (x, y) in sorted(groups, key=lambda k: (-k[1], k[0])):
        names = ' = '.join(groups[(x, y)])
        left = x > 0.55
        below = any(abs(x - px) < 0.10 and abs(y - py) < 0.030
                    for px, py in placed)
        ax.annotate(names, (x, y), fontsize=7,
                    ha='right' if left else 'left',
                    va='top' if below else 'bottom',
                    xytext=(-4 if left else 4, -4 if below else 4),
                    textcoords='offset points')
        placed.append((x, y))
    ax.plot([0, 1], [0, 1], ls=':', c='#999999', lw=1)
    ax.axhline(0.90, ls='--', c='#bbbbbb', lw=0.8)
    ax.axvline(0.90, ls='--', c='#bbbbbb', lw=0.8)
    ax.set_xlabel('separability, oil vs clean soil')
    ax.set_ylabel('separability, oil vs burn scar')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title('No index does both jobs')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_FIGURES, 'fig_dissociation.png'), dpi=200)
    plt.close(fig)

    # Figure 2: per-event spread for the named features. Any feature
    # missing from the data is dropped with a warning rather than
    # raising, so a reduced run still produces a figure.
    show = [f for f in FIG2_DISCRIMINATORS + FIG2_DETECTORS if f in ev_cols]
    missing = [f for f in FIG2_DISCRIMINATORS + FIG2_DETECTORS
               if f not in ev_cols]
    if missing:
        print('Figure 2: not in data, skipped: %s' % ', '.join(missing))
    n_disc = sum(1 for f in FIG2_DISCRIMINATORS if f in ev_cols)
    data = []
    for f in show:
        vals = [separability(oil_c[f], a) for a in ev_cols[f].values()]
        data.append([v for v in vals if np.isfinite(v)])
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.boxplot(data, tick_labels=[pretty(f) for f in show], showfliers=False)
    for i, d in enumerate(data, start=1):
        ax.scatter(np.full(len(d), i) + rng.uniform(-0.12, 0.12, len(d)),
                   d, s=9, alpha=0.45, c='#666666', edgecolor='none')
    ax.axhline(0.80, ls='--', c='#bbbbbb', lw=0.8)
    # Divider between the two groups, so the reader does not have to
    # take the grouping on trust from the caption alone.
    if 0 < n_disc < len(show):
        ax.axvline(n_disc + 0.5, ls='-', c='#333333', lw=0.9)
        ax.text(n_disc / 2.0 + 0.5, 1.02, 'discriminate',
                ha='center', va='bottom', fontsize=8, transform=
                ax.get_xaxis_transform())
        ax.text((n_disc + len(show)) / 2.0 + 0.5, 1.02, 'detect',
                ha='center', va='bottom', fontsize=8, transform=
                ax.get_xaxis_transform())
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel('separability, oil vs burn')
    # Padded so the group labels above the axis do not meet the title.
    ax.set_title('Per burn event (n = %d)' % len(big_events), pad=20)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_FIGURES, 'fig_per_event.png'), dpi=200)
    plt.close(fig)

    # Figure 3: the two-dimensional feature space
    fig, ax = plt.subplots(figsize=(6.6, 6))
    for rs, lab, col, a, z in [
            (burn, 'burn', '#c0392b', 0.18, 1),
            (clean, 'clean soil', '#2980b9', 0.25, 2),
            (oil, 'oil', '#f1c40f', 0.95, 3)]:
        ax.scatter(column(rs, det), column(rs, dis), s=10, alpha=a,
                   c=col, label=lab, edgecolor='none', zorder=z)

    # The fitted discriminant, drawn so that the 0.960 of the text has a
    # visible counterpart. This is the full-sample fit, not one of the
    # cross-validation folds: the folds differ only in which fire is
    # absent, and drawing 260 nearly coincident lines would say less.
    ob = np.column_stack([column(oil, det), column(oil, dis)])
    bb = np.column_stack([column(burn, det), column(burn, dis)])
    ob = ob[clean_rows([ob[:, 0], ob[:, 1]])]
    bb = bb[clean_rows([bb[:, 0], bb[:, 1]])]
    if len(ob) >= MIN_N and len(bb) >= MIN_N:
        wv = fisher_ld(bb, ob)
        # Threshold midway between the projected class means.
        thr = 0.5 * (ob.mean(axis=0) @ wv + bb.mean(axis=0) @ wv)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        if abs(wv[1]) > abs(wv[0]):
            xs = np.array([x0, x1])
            ys = (thr - wv[0] * xs) / wv[1]
        else:
            ys = np.array([y0, y1])
            xs = (thr - wv[1] * ys) / wv[0]
        ax.plot(xs, ys, color='#222222', lw=1.4, ls='--', zorder=4,
                label='oil / burn boundary')
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)

    ax.set_xlabel(pretty(det))
    ax.set_ylabel(pretty(dis))
    ax.legend(frameon=False)
    ax.set_title('%s separates oil from soil; %s from char'
                 % (pretty(det), pretty(dis)))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_FIGURES, 'fig_feature_space.png'), dpi=200)
    plt.close(fig)

    print('written:')
    print('  %s/separability.csv' % OUT_RESULTS)
    print('  %s/fig_dissociation.png' % OUT_FIGURES)
    print('  %s/fig_per_event.png' % OUT_FIGURES)
    print('  %s/fig_feature_space.png' % OUT_FIGURES)
    print()
    print('CAVEAT, to be carried into the manuscript: the oil class is')
    print('one plume, one site, one date. Every interval quoted for the')
    print('oil side is a within-plume interval and says nothing about')
    print('how these indices behave on other spills.')


if __name__ == '__main__':
    main()
