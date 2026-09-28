#!/usr/bin/env python3
"""
Recover the plume polygon from the spectra.

The hand-digitised polygon was never saved and the sample export dropped
pixel geometries, so the label exists nowhere. It can be rebuilt: the
ten-band reflectance tuple is unique across every pixel sampled at the
Narlı site, so matching those tuples against a full export of the 3 x 3
km window relocates each sampled pixel exactly.

Inputs
    oil_window_pixels.csv   every pixel in OIL_WINDOW, with coordinates,
                            exported by the RECOVER block appended to
                            narli_burn_confusion_gee.js
    oil_burn_samples.csv    the data of record

Outputs
    oil_plume_polygon.geojson   dissolved outline of the oil class
    oil_plume_pixels.csv        matched pixel centres, one row per pixel

The dissolve is exact rather than traced: each matched pixel contributes
its 20 m cell, shared edges between adjacent cells cancel, and the
surviving edges are chained into rings. The result is the true footprint
of the sampled class, not a smoothed approximation of it.

Usage:
    python3 recover_plume.py [oil_window_pixels.csv] [oil_burn_samples.csv]
"""

import csv
import json
import os
import sys
from collections import Counter, defaultdict

import utm

BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']

# The archive stores reflectance to four decimals. The recovery export
# may carry more, so both sides are rounded to the archive's precision
# before matching.
DP = 4

# Classes to recover, and which one defines the plume polygon.
CLASSES = ['oil', 'oil_edge', 'control_clean']
PLUME_CLASS = 'oil'

PIXEL_M = 20.0


def key(row):
    try:
        return tuple(round(float(row[b]), DP) for b in BANDS)
    except (TypeError, ValueError):
        return None


def rings_from_cells(cells):
    """
    Dissolve a set of (col, row) grid cells into closed rings.

    Each cell contributes four directed edges counter-clockwise. An edge
    shared by two adjacent cells appears once in each direction and the
    pair cancels, so only the outer boundary survives. The survivors are
    then chained head to tail into rings.
    """
    edges = set()
    for (i, j) in cells:
        for a, b in (((i, j), (i + 1, j)),
                     ((i + 1, j), (i + 1, j + 1)),
                     ((i + 1, j + 1), (i, j + 1)),
                     ((i, j + 1), (i, j))):
            if (b, a) in edges:
                edges.remove((b, a))
            else:
                edges.add((a, b))

    nxt = defaultdict(list)
    for a, b in edges:
        nxt[a].append(b)

    rings = []
    remaining = set(edges)
    while remaining:
        a, b = next(iter(remaining))
        ring = [a]
        cur_from, cur_to = a, b
        while True:
            remaining.discard((cur_from, cur_to))
            ring.append(cur_to)
            if cur_to == a:
                break
            outs = [t for t in nxt[cur_to] if (cur_to, t) in remaining]
            if not outs:
                break
            cur_from, cur_to = cur_to, outs[0]
        rings.append(ring)
    return rings


def signed_area(ring):
    s = 0.0
    for (x1, y1), (x2, y2) in zip(ring, ring[1:]):
        s += x1 * y2 - x2 * y1
    return s / 2.0


def point_in_ring(pt, ring):
    x, y = pt
    inside = False
    for (x1, y1), (x2, y2) in zip(ring, ring[1:]):
        if (y1 > y) != (y2 > y):
            xin = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x < xin:
                inside = not inside
    return inside


def main():
    win_path = sys.argv[1] if len(sys.argv) > 1 else 'oil_window_pixels.csv'
    smp_path = sys.argv[2] if len(sys.argv) > 2 else 'oil_burn_samples.csv'
    for p in (win_path, smp_path):
        if not os.path.exists(p):
            raise SystemExit(
                'Not found: %s\n'
                'Export the window first: append recover_plume_gee.js to '
                'narli_burn_confusion_gee.js, set RUN_STAGE = \'RECOVER\', '
                'run it, and download the CSV here.' % p)

    with open(win_path, newline='') as fh:
        win = list(csv.DictReader(fh))
    with open(smp_path, newline='') as fh:
        smp = [r for r in csv.DictReader(fh)
               if r['class_name'] in CLASSES]

    print('window pixels exported : %d' % len(win))
    print('site pixels to relocate: %d' % len(smp))

    # Collisions would make a match ambiguous, so they are excluded
    # rather than resolved arbitrarily.
    counts = Counter(key(r) for r in win)
    ambiguous = {k for k, c in counts.items() if c > 1 and k is not None}
    lookup = {}
    for r in win:
        k = key(r)
        if k is not None and k not in ambiguous:
            lookup[k] = (float(r['longitude']), float(r['latitude']))
    print('ambiguous band tuples  : %d (excluded)' % len(ambiguous))

    matched, missed = [], []
    for r in smp:
        k = key(r)
        if k in lookup:
            lon, lat = lookup[k]
            matched.append((r['class_name'], lon, lat))
        else:
            missed.append(r['class_name'])
    print()
    print('matched : %d of %d (%.1f %%)'
          % (len(matched), len(smp), 100.0 * len(matched) / max(len(smp), 1)))
    if missed:
        print('unmatched by class: %s' % dict(Counter(missed)))
        print('  A shortfall means the recovery export did not reproduce '
              'the archived composite. Check that the printed pre_date and '
              'post_date match the archive before trusting the result.')

    with open('oil_plume_pixels.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['class_name', 'longitude', 'latitude'])
        w.writerows(matched)

    # --- grid, in UTM ------------------------------------------------
    # The Sentinel-2 sampling grid is UTM, not geographic. At this
    # latitude meridian convergence rotates it about 1.2 degrees from
    # true north, so in longitude and latitude no two of the 22380
    # window pixels share a coordinate and no grid can be recovered
    # there. Projected to the image's own UTM zone it is exactly regular
    # at 20 m, which is what makes the dissolve below meaningful.
    zone = utm.zone_from_lon(
        sorted(float(r['longitude']) for r in win)[len(win) // 2])
    xs, ys = [], []
    for r in win:
        x, y = utm.forward(float(r['longitude']), float(r['latitude']), zone)
        xs.append(x)
        ys.append(y)
    x0, y0 = min(xs), min(ys)

    def cell_of(lon, lat):
        x, y = utm.forward(lon, lat, zone)
        return (int(round((x - x0) / PIXEL_M)), int(round((y - y0) / PIXEL_M)))

    def corner(i, j):
        return utm.inverse(x0 + (i - 0.5) * PIXEL_M,
                           y0 + (j - 0.5) * PIXEL_M, zone)

    # The grid must account for every exported pixel exactly once. If it
    # does not, the projection or the pixel size is wrong and every
    # number below would be quietly meaningless.
    occupied = {cell_of(float(r['longitude']), float(r['latitude']))
                for r in win}
    print()
    print('UTM zone    : %d' % zone)
    print('grid cells  : %d occupied by %d exported pixels' % (len(occupied), len(win)))
    if len(occupied) != len(win):
        raise SystemExit(
            'Grid collision: %d pixels fell into %d cells. The sampling '
            'grid is not %g m in UTM zone %d, so the polygon would be '
            'wrong. Stopping rather than writing it.'
            % (len(win), len(occupied), PIXEL_M, zone))
    span_i = max(i for i, _ in occupied) - min(i for i, _ in occupied) + 1
    span_j = max(j for _, j in occupied) - min(j for _, j in occupied) + 1
    print('grid extent : %d x %d cells (%.1f x %.1f km)'
          % (span_i, span_j, span_i * PIXEL_M / 1000, span_j * PIXEL_M / 1000))

    cells = {cell_of(lon, lat) for c, lon, lat in matched if c == PLUME_CLASS}
    print()
    print('plume cells : %d' % len(cells))
    print('area        : %.2f ha  (%d cells x %.0f m x %.0f m)'
          % (len(cells) * PIXEL_M * PIXEL_M / 1e4, len(cells), PIXEL_M, PIXEL_M))

    rings = [[corner(i, j) for (i, j) in ring]
             for ring in rings_from_cells(cells)]
    outers = [r for r in rings if signed_area(r) > 0]
    holes = [r for r in rings if signed_area(r) <= 0]
    print('rings       : %d outer, %d holes' % (len(outers), len(holes)))

    polys = []
    for o in outers:
        mine = [h for h in holes if point_in_ring(h[0], o)]
        polys.append([o] + mine)

    gj = {
        'type': 'Feature',
        'properties': {
            'name': 'Narli oil plume, recovered from sampled spectra',
            'class': PLUME_CLASS,
            'n_pixels': len(cells),
            'area_ha': round(len(cells) * PIXEL_M * PIXEL_M / 1e4, 4),
            'pixel_size_m': PIXEL_M,
            'note': ('Reconstructed by matching the ten-band reflectance '
                     'tuple of each sampled pixel against a full export of '
                     'the 3 x 3 km analysis window. This is the exact '
                     'footprint of the sampled oil class, not the original '
                     'hand-drawn outline, which was not archived.')
        },
        'geometry': {'type': 'MultiPolygon', 'coordinates': polys}
        if len(polys) != 1 else {'type': 'Polygon', 'coordinates': polys[0]}
    }
    with open('oil_plume_polygon.geojson', 'w') as fh:
        json.dump(gj, fh, indent=1)

    print()
    print('written: oil_plume_polygon.geojson, oil_plume_pixels.csv')


if __name__ == '__main__':
    main()
