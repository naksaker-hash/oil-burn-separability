#!/usr/bin/env python3
"""
Nested leave-one-burn-event-out validation of the two-feature
discriminant.

Why this exists
---------------
In oil_burn_separability.py the detector and the discriminator are
chosen from the ranking computed on the complete sample, and only then
is the pair validated leave-one-event-out. The held-out fire therefore
contributed to the choice of features, which makes the reported figure
optimistic by an unknown amount. A reviewer asked, correctly, whether
the selection happens inside the fold.

Here it does. For every held-out fire:

    1. remove that fire entirely
    2. using only the remaining fires and the training half of the oil
       pixels, pick the strongest detector (oil vs clean soil) and the
       strongest discriminator (oil vs burn)
    3. fit Fisher's linear discriminant on those two features
    4. score the held-out fire against the withheld oil pixels

Baselines are evaluated on exactly the same folds so that the pair can
be compared with each index alone, and with the fixed pair used in the
manuscript, on equal terms.

The oil side cannot be held out: there is one spill. Oil pixels are
split at random, which does not confer independence on them. What the
folds estimate is generalisation to UNSEEN FIRES.

Dependencies: numpy only.

Usage:
    python3 nested_loeo.py [path/to/oil_burn_samples.csv]
"""

import csv
import os
import sys
import numpy as np

from oil_burn_separability import (
    load, subset, column, separability, fisher_ld, clean_rows,
    MIN_N, MIN_EVENT_N, BOOT_SEED, LOEO_REPEATS, OUT_RESULTS, pretty,
)


def pctile(a, q):
    return np.percentile(a, q) if a.size else np.nan


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'oil_burn_samples.csv'
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)
    os.makedirs(OUT_RESULTS, exist_ok=True)

    rows, feats = load(path)
    oil = subset(rows, 'oil')
    clean = subset(rows, 'control_clean')
    burn = subset(rows, 'burn')

    big = {}
    for e in sorted({r['event_id'] for r in burn}):
        b = [r for r in burn if r['event_id'] == e]
        if len(b) >= MIN_EVENT_N:
            big[e] = b

    oil_c = {f: column(oil, f) for f in feats}
    clean_c = {f: column(clean, f) for f in feats}
    ev_c = {f: {k: column(v, f) for k, v in big.items()} for f in feats}

    # Features usable at all: finite in the oil class and in clean soil.
    cand = [f for f in feats
            if np.isfinite(oil_c[f]).sum() >= MIN_N
            and np.isfinite(clean_c[f]).sum() >= MIN_N]

    keys = list(big)
    print('=' * 70)
    print('NESTED LEAVE-ONE-BURN-EVENT-OUT VALIDATION')
    print('=' * 70)
    print('source            : %s' % path)
    print('candidate features: %d' % len(cand))
    print('burn events       : %d (>= %d pixels)' % (len(keys), MIN_EVENT_N))
    print('oil pixels        : %d' % len(oil))
    print('oil-split repeats : %d' % LOEO_REPEATS)
    print()

    rng = np.random.default_rng(BOOT_SEED)
    n_oil = len(oil)

    scores = {k: [] for k in
              ('nested_pair', 'nested_det', 'nested_dis',
               'fixed_pair', 'fixed_det', 'fixed_dis')}
    picked_det, picked_dis = {}, {}
    per_fold = []

    FIXED_DET, FIXED_DIS = 'd_NBR', 'MIRBI'

    for rep in range(LOEO_REPEATS):
        idx = rng.permutation(n_oil)
        half = n_oil // 2
        tr_i, te_i = idx[:half], idx[half:]

        # Detector selection uses oil-train against clean soil. No burn
        # event enters it, so it cannot leak the held-out fire; it is
        # recomputed per repeat anyway because the oil half changes.
        det_rank = []
        for f in cand:
            s = separability(oil_c[f][tr_i], clean_c[f])
            if np.isfinite(s):
                det_rank.append((s, f))
        det_rank.sort(reverse=True)
        det = det_rank[0][1]

        for k in keys:
            # --- inner selection, held-out fire excluded -------------
            dis, best = None, -1.0
            for f in cand:
                pooled = np.concatenate([ev_c[f][kk] for kk in keys if kk != k])
                s = separability(oil_c[f][tr_i], pooled)
                if np.isfinite(s) and s > best:
                    best, dis = s, f
            if dis is None:
                continue
            picked_det[det] = picked_det.get(det, 0) + 1
            picked_dis[dis] = picked_dis.get(dis, 0) + 1

            row = {'rep': rep, 'event_id': k, 'det': det, 'dis': dis}

            for tag, a, b in (('nested', det, dis),
                              ('fixed', FIXED_DET, FIXED_DIS)):
                if a not in oil_c or b not in oil_c:
                    continue
                o2 = np.column_stack([oil_c[a], oil_c[b]])
                m = clean_rows([o2[:, 0], o2[:, 1]])
                otr = o2[np.intersect1d(tr_i, np.flatnonzero(m))]
                ote = o2[np.intersect1d(te_i, np.flatnonzero(m))]

                def ev2(kk):
                    x = np.column_stack([ev_c[a][kk], ev_c[b][kk]])
                    return x[clean_rows([x[:, 0], x[:, 1]])]

                tr = [ev2(kk) for kk in keys if kk != k]
                tr = [t for t in tr if len(t) >= MIN_N]
                bte = ev2(k)
                if not tr or len(bte) < MIN_N or len(otr) < MIN_N \
                        or len(ote) < MIN_N:
                    continue
                btr = np.vstack(tr)
                w = fisher_ld(btr, otr)
                s_pair = separability(ote @ w, bte @ w)
                s_det = separability(ote[:, 0], bte[:, 0])
                s_dis = separability(ote[:, 1], bte[:, 1])
                for name, val in (('%s_pair' % tag, s_pair),
                                  ('%s_det' % tag, s_det),
                                  ('%s_dis' % tag, s_dis)):
                    if np.isfinite(val):
                        scores[name].append(val)
                        row[name] = val
            per_fold.append(row)

    print('SELECTED INSIDE THE FOLD')
    print('  detector      : %s' %
          ', '.join('%s %d/%d' % (pretty(f), n, len(per_fold))
                    for f, n in sorted(picked_det.items(), key=lambda t: -t[1])))
    print('  discriminator : %s' %
          ', '.join('%s %d/%d' % (pretty(f), n, len(per_fold))
                    for f, n in sorted(picked_dis.items(), key=lambda t: -t[1])))
    print()

    print('SEPARABILITY ON WITHHELD FIRES')
    print('%-32s %6s %7s %7s %7s' % ('model', 'folds', 'median', 'p10', 'p90'))
    print('-' * 64)
    label = {
        'nested_pair': 'nested pair (selected in fold)',
        'nested_det':  '  its detector alone',
        'nested_dis':  '  its discriminator alone',
        'fixed_pair':  'fixed pair (%s + %s)' % (pretty(FIXED_DET), pretty(FIXED_DIS)),
        'fixed_det':   '  %s alone' % pretty(FIXED_DET),
        'fixed_dis':   '  %s alone' % pretty(FIXED_DIS),
    }
    for k in ('nested_pair', 'nested_det', 'nested_dis',
              'fixed_pair', 'fixed_det', 'fixed_dis'):
        a = np.array(scores[k])
        print('%-32s %6d %7.3f %7.3f %7.3f' %
              (label[k], a.size, np.median(a) if a.size else np.nan,
               pctile(a, 10), pctile(a, 90)))
    print()

    npair = np.array(scores['nested_pair'])
    fpair = np.array(scores['fixed_pair'])
    if npair.size and fpair.size:
        print('nested median minus fixed median : %+.4f'
              % (np.median(npair) - np.median(fpair)))
        print('fire-level test units            : %d' % len(keys))
        print('evaluations (units x repeats)    : %d' % npair.size)
    print()

    out = os.path.join(OUT_RESULTS, 'nested_loeo.csv')
    cols = ['rep', 'event_id', 'det', 'dis', 'nested_pair', 'nested_det',
            'nested_dis', 'fixed_pair', 'fixed_det', 'fixed_dis']
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerows(per_fold)
    print('written: %s' % out)


if __name__ == '__main__':
    main()
