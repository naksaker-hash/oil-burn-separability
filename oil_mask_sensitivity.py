#!/usr/bin/env python3
"""
Sensitivity of the main result to how the oil plume is delineated.

The oil reference is a hand-drawn polygon on a SWIR composite. A
reviewer asked, reasonably, whether the separability depends on where
that line was put. This repeats the univariate analysis under three
masks of increasing inclusiveness:

    core      the sampled oil pixels whose four 20 m neighbours are all
              themselves oil -- the interior, with the digitised edge
              stripped off
    plume     the 218 oil pixels as sampled, i.e. the manuscript's mask
    inclusive plume plus the 100 oil_edge pixels, i.e. the boundary read
              generously outwards

The three masks are nested, so the comparison isolates the effect of
the boundary rather than of a different sample.

Pixel positions come from oil_plume_pixels.csv, which recover_plume.py
writes in the same row order as the class-filtered sample table, so the
two join by position. Positions are projected to UTM before neighbours
are counted: the Sentinel-2 grid is UTM, and at this latitude it is
rotated about a degree from true north, so a neighbour test in degrees
finds nothing.

Dependencies: numpy only.

Usage:
    python3 oil_mask_sensitivity.py
"""

import csv
import os
import sys
import numpy as np

from oil_burn_separability import (
    load, subset, column, separability, MIN_N, MIN_EVENT_N, OUT_RESULTS,
    boot_ci_events, pretty,
)
from utm import forward as utm_forward

PIX = 20.0          # Sentinel-2 SWIR sampling, metres
TOL = 4.0           # neighbour match tolerance, metres
REPORT = ['d_NBR', 'd_FMR', 'NBR', 'd_logBAI', 'd_BAI',
          'MIRBI', 'NBR2', 'NDVI', 'B11', 'B12']


def grid_cells(lons, lats):
    """(col, row) indices on the local 20 m UTM grid."""
    xs, ys = [], []
    zone = int((np.median(lons) + 180.0) // 6) + 1
    for lon, lat in zip(lons, lats):
        x, y = utm_forward(lon, lat, zone)
        xs.append(x)
        ys.append(y)
    xs, ys = np.array(xs), np.array(ys)
    col = np.round((xs - xs.min()) / PIX).astype(int)
    row = np.round((ys - ys.min()) / PIX).astype(int)
    err = max(np.abs((xs - xs.min()) / PIX - col).max(),
              np.abs((ys - ys.min()) / PIX - row).max()) * PIX
    if err > TOL:
        raise SystemExit('pixel centres do not lie on a %.0f m grid '
                         '(worst residual %.1f m)' % (PIX, err))
    return list(zip(col, row))


def main():
    smp_path = 'oil_burn_samples.csv'
    pix_path = 'oil_plume_pixels.csv'
    for p in (smp_path, pix_path):
        if not os.path.exists(p):
            raise SystemExit('Not found: %s' % p)

    rows, feats = load(smp_path)
    with open(pix_path, newline='') as fh:
        pix = list(csv.DictReader(fh))

    order = [r for r in rows if r['class_name'] in
             ('oil', 'oil_edge', 'control_clean')]
    if len(order) != len(pix):
        raise SystemExit('%d sample rows but %d located pixels; the join is '
                         'positional and cannot be trusted'
                         % (len(order), len(pix)))
    for a, b in zip(order, pix):
        if a['class_name'] != b['class_name']:
            raise SystemExit('class mismatch in positional join')

    lon = np.array([float(r['longitude']) for r in pix])
    lat = np.array([float(r['latitude']) for r in pix])
    cells = grid_cells(lon, lat)

    is_oil = np.array([r['class_name'] == 'oil' for r in order])
    is_edge = np.array([r['class_name'] == 'oil_edge' for r in order])

    oil_cells = {c for c, k in zip(cells, is_oil) if k}
    core_mask = np.array([
        k and all((c[0] + dx, c[1] + dy) in oil_cells
                  for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)))
        for c, k in zip(cells, is_oil)])

    oil_rows = [r for r in order if r['class_name'] == 'oil']
    edge_rows = [r for r in order if r['class_name'] == 'oil_edge']
    core_rows = [r for r, k in zip(order, core_mask) if k]

    masks = [
        ('core', core_rows),
        ('plume', oil_rows),
        ('inclusive', oil_rows + edge_rows),
    ]

    clean = subset(rows, 'control_clean')
    burn = subset(rows, 'burn')
    big = {}
    for e in sorted({r['event_id'] for r in burn}):
        b = [r for r in burn if r['event_id'] == e]
        if len(b) >= MIN_EVENT_N:
            big[e] = b

    print('=' * 70)
    print('OIL-MASK SENSITIVITY')
    print('=' * 70)
    for name, rs in masks:
        print('  %-10s %4d pixels  (%.2f ha)' % (name, len(rs),
                                                 len(rs) * PIX * PIX / 1e4))
    print()
    print('%-12s %s' % ('feature',
          ''.join('%22s' % ('%s vs clean / burn' % n) for n, _ in masks)))
    print('-' * 80)

    out = []
    for f in REPORT:
        if f not in feats:
            continue
        line = '%-12s' % pretty(f)
        rec = {'feature': f}
        for name, rs in masks:
            o = column(rs, f)
            sc = separability(o, column(clean, f))
            sb = separability(o, column(burn, f))
            rec['%s_clean' % name] = sc
            rec['%s_burn' % name] = sb
            line += '%22s' % ('%.3f / %.3f' % (sc, sb))
        out.append(rec)
        print(line)
    print()

    # Interval on the headline discriminator under each mask, so that
    # the masks can be compared against the width of the uncertainty
    # rather than against zero.
    ev_c = {}
    for f in ('MIRBI', 'NBR2', 'd_logBAI'):
        if f not in feats:
            continue
        ev_c[f] = {k: column(v, f) for k, v in big.items()}
    print('CLUSTER-BOOTSTRAP INTERVAL (burn events resampled)')
    print('%-12s %-10s %8s   %s' % ('feature', 'mask', 'S', '95% CI'))
    print('-' * 52)
    for f in ev_c:
        for name, rs in masks:
            lo, hi = boot_ci_events(column(rs, f), ev_c[f])
            print('%-12s %-10s %8.3f   [%.3f, %.3f]' %
                  (pretty(f), name, separability(column(rs, f),
                                                 column(burn, f)), lo, hi))
    print()

    os.makedirs(OUT_RESULTS, exist_ok=True)
    p = os.path.join(OUT_RESULTS, 'oil_mask_sensitivity.csv')
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print('written: %s' % p)


if __name__ == '__main__':
    main()
