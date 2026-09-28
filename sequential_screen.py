#!/usr/bin/env python3
"""
A two-stage screen: detect first, then reject fire.

The separability analysis measures how far apart two classes sit along
a feature. It does not set a threshold, and it does not say what a
detector would do when it meets a landscape containing clean soil,
unburned vegetation and fire scars at once. This does, in the simplest
form that can be written down.

    Stage 1  flag candidates on the detector, at a threshold set to
             retain 90 per cent of training oil pixels against the
             background classes
    Stage 2  reject candidates on the discriminator, at a threshold set
             to retain 90 per cent of the stage 1 survivors

Both thresholds are fitted on training data only; the fire under test
is withheld entirely, together with its own control annulus and half
the oil pixels. Thresholds are quantiles of the oil training half, so
they are set by the target class rather than tuned against the outcome.

Reported per fold: what fraction of withheld oil, withheld burn and
withheld background survives each stage. Precision is then computed at
an assumed scene composition, which is an illustration and not a
measurement -- no scene is classified here.

Dependencies: numpy only.

Usage:
    python3 sequential_screen.py [path/to/oil_burn_samples.csv]
"""

import csv
import os
import sys
import numpy as np

from oil_burn_separability import (
    load, subset, column, clean_rows, MIN_N, MIN_EVENT_N, BOOT_SEED,
    LOEO_REPEATS, OUT_RESULTS, pretty,
)

DET = 'd_NBR'
DIS = 'MIRBI'
KEEP = 0.90          # fraction of training oil retained at each stage
PREVALENCE = 0.0097  # 8.72 ha plume in the 3 x 3 km analysis window


def side(a, b):
    """+1 if class a sits above class b on this feature, else -1."""
    return 1.0 if np.nanmedian(a) > np.nanmedian(b) else -1.0


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'oil_burn_samples.csv'
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)
    rows, feats = load(path)
    oil = subset(rows, 'oil')
    clean = subset(rows, 'control_clean')
    burn = subset(rows, 'burn')
    ctrl = subset(rows, 'control_unburned')

    ev = sorted({r['event_id'] for r in burn})
    big = {}
    for e in ev:
        b = [r for r in burn if r['event_id'] == e]
        if len(b) >= MIN_EVENT_N:
            big[e] = b
    keys = list(big)

    def two(rs):
        m = np.column_stack([column(rs, DET), column(rs, DIS)])
        return m[clean_rows([m[:, 0], m[:, 1]])]

    O = two(oil)
    C = two(clean)
    B = {k: two(v) for k, v in big.items()}
    K = {}
    for k in keys:
        rs = [r for r in ctrl if r['event_id'] == k]
        if rs:
            K[k] = two(rs)

    # Orientation is fixed once, on the pooled training-side medians:
    # which way is "oil-like" is a property of the feature, not of the
    # fold, and re-deciding it per fold would let a fold flip the rule.
    bg_all = np.vstack([C] + [v for v in K.values()])
    s_det = side(O[:, 0], bg_all[:, 0])
    s_dis = side(O[:, 1], np.vstack(list(B.values()))[:, 1])

    print('=' * 70)
    print('SEQUENTIAL SCREEN: detect, then reject fire')
    print('=' * 70)
    print('stage 1 feature : %s  (oil %s background)'
          % (pretty(DET), 'above' if s_det > 0 else 'below'))
    print('stage 2 feature : %s  (oil %s burn)'
          % (pretty(DIS), 'above' if s_dis > 0 else 'below'))
    print('retention       : %.0f%% of training oil at each stage'
          % (100 * KEEP))
    print('test units      : %d fires x %d oil splits' % (len(keys),
                                                          LOEO_REPEATS))
    print()

    rng = np.random.default_rng(BOOT_SEED)
    acc = {k: [] for k in ('oil1', 'oil2', 'burn1', 'burn2',
                           'clean1', 'clean2', 'ctrl1', 'ctrl2')}
    per_fold = []

    for rep in range(LOEO_REPEATS):
        idx = rng.permutation(len(O))
        half = len(O) // 2
        Otr, Ote = O[idx[:half]], O[idx[half:]]
        # Thresholds depend only on the oil training half, so they are
        # constant within a repetition.
        t1 = np.quantile(s_det * Otr[:, 0], 1.0 - KEEP)
        surv = Otr[s_det * Otr[:, 0] >= t1]
        if len(surv) < MIN_N:
            continue
        t2 = np.quantile(s_dis * surv[:, 1], 1.0 - KEEP)

        def rate(X):
            if len(X) == 0:
                return np.nan, np.nan
            p1 = s_det * X[:, 0] >= t1
            p2 = p1 & (s_dis * X[:, 1] >= t2)
            return p1.mean(), p2.mean()

        for k in keys:
            bg_clean = C
            bg_ctrl = K.get(k, np.empty((0, 2)))
            o1, o2 = rate(Ote)
            b1, b2 = rate(B[k])
            c1, c2 = rate(bg_clean)
            k1, k2 = rate(bg_ctrl)
            for name, v in (('oil1', o1), ('oil2', o2), ('burn1', b1),
                            ('burn2', b2), ('clean1', c1), ('clean2', c2),
                            ('ctrl1', k1), ('ctrl2', k2)):
                if np.isfinite(v):
                    acc[name].append(v)
            per_fold.append({'rep': rep, 'event_id': k, 'oil1': o1,
                             'oil2': o2, 'burn1': b1, 'burn2': b2,
                             'clean1': c1, 'clean2': c2, 'ctrl1': k1,
                             'ctrl2': k2})

    def med(k):
        a = np.array(acc[k])
        return np.median(a) if a.size else np.nan

    print('FRACTION SURVIVING, MEDIAN OVER FOLDS')
    print('%-24s %12s %12s' % ('class', 'after st.1', 'after st.2'))
    print('-' * 50)
    for lab, a, b in (('oil (withheld)', 'oil1', 'oil2'),
                      ('burn (withheld fire)', 'burn1', 'burn2'),
                      ('clean soil', 'clean1', 'clean2'),
                      ('unburned control', 'ctrl1', 'ctrl2')):
        print('%-24s %12.3f %12.3f' % (lab, med(a), med(b)))
    print()

    # Illustrative precision. Scene composition is assumed, not measured.
    tpr1, tpr2 = med('oil1'), med('oil2')
    bg1 = np.nanmean([med('clean1'), med('ctrl1')])
    bg2 = np.nanmean([med('clean2'), med('ctrl2')])
    print('ILLUSTRATIVE PRECISION AT P = %.4f' % PREVALENCE)
    print('%-28s %10s %10s' % ('burnt fraction of scene', 'one stage',
                               'two stages'))
    print('-' * 52)
    for burn_frac in (0.01, 0.05, 0.10):
        rest = 1.0 - PREVALENCE - burn_frac
        fp1 = med('burn1') * burn_frac + bg1 * rest
        fp2 = med('burn2') * burn_frac + bg2 * rest
        p1 = tpr1 * PREVALENCE / (tpr1 * PREVALENCE + fp1)
        p2 = tpr2 * PREVALENCE / (tpr2 * PREVALENCE + fp2)
        print('%-28s %10.3f %10.3f' % ('%.0f%%' % (100 * burn_frac), p1, p2))
    print()
    print('Composition is assumed. No scene is classified here, and these')
    print('are not measured accuracies.')
    print()

    os.makedirs(OUT_RESULTS, exist_ok=True)
    p = os.path.join(OUT_RESULTS, 'sequential_screen.csv')
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_fold[0].keys()))
        w.writeheader()
        w.writerows(per_fold)
    print('written: %s' % p)


if __name__ == '__main__':
    main()
