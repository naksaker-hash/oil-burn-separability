#!/usr/bin/env python3
"""
Which two burn events Figure 4 shows, and where they are.

The Earth Engine script picks the panels by rule, not by eye: the two
events nearest the plume in normalised (pre-fire NDVI, slope) space,
each axis divided by a nominal scale of 0.1 NDVI units and 2 degrees.
Choosing by appearance would let the figure argue something the data
does not.

The rule is reproduced here from fig05_event_points.csv so the caption
can be written without re-running Earth Engine, and so the choice can be
shown to be insensitive to which plume slope value is used -- the one
disputed between this study and paper 1.

Also prints the distances behind section 4.5's claim that substrate
matching selected climatically comparable terrain without being told to.

Usage:
    python3 fig04_panel_choice.py [fig05_event_points.csv]
"""

import csv
import math
import os
import sys

import numpy as np

from provinces import repair

# Rupture point, from make_figs_45.py.
SPILL_LON, SPILL_LAT = 37.130123, 37.350183

# Nominal scales from make_figs_45_gee.js. The two axes have different
# units and ranges, so each difference is divided before combining.
NDVI_SCALE = 0.1
SLOPE_SCALE = 2.0

# Plume pre-fire NDVI, agreed by every route.
PLUME_NDVI = 0.193

# Plume slope is not agreed: the per-pixel median over the sampled plume
# gives 1.00 deg, the Earth Engine 200 m buffer mean gives 0.799, and
# paper 1 reports 0.69. All three are tried, because a panel choice that
# depended on which one is used would not be worth defending.
SLOPE_CANDIDATES = [
    (1.00, 'per-pixel median over the sampled plume'),
    (0.799, 'Earth Engine 200 m buffer mean'),
    (0.69, 'as paper 1 reports'),
]


def haversine_km(lon, lat):
    R = 6371.0088
    p1, p2 = math.radians(SPILL_LAT), math.radians(lat)
    dp = p2 - p1
    dl = math.radians(lon - SPILL_LON)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


def ranked(rows, plume_ndvi, plume_slope):
    out = []
    for r in rows:
        dn = (float(r['pre_ndvi']) - plume_ndvi) / NDVI_SCALE
        ds = (float(r['slope_mean']) - plume_slope) / SLOPE_SCALE
        out.append((math.hypot(dn, ds), r))
    out.sort(key=lambda t: t[0])
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'fig05_event_points.csv'
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))

    print('PANEL SELECTION, two nearest in normalised (NDVI, slope) space')
    print('events considered: %d' % len(rows))
    picks = []
    for slope, note in SLOPE_CANDIDATES:
        top = ranked(rows, PLUME_NDVI, slope)[:2]
        picks.append([r['event_id'] for _, r in top])
        print('\n  plume slope %.3f deg (%s)' % (slope, note))
        for i, (d, r) in enumerate(top, start=1):
            print('    burn %d  %-30s dist %.3f' % (i, r['event_id'], d))

    if all(p == picks[0] for p in picks):
        print('\n  Same two events, same order, under every slope value. '
              'The choice does not depend on the disputed number.')
    else:
        print('\n  WARNING: the panel choice changes with the slope value. '
              'Resolve the slope before writing the caption.')

    print('\nPANEL DETAIL')
    top = ranked(rows, PLUME_NDVI, SLOPE_CANDIDATES[0][0])[:2]
    for i, (d, r) in enumerate(top, start=1):
        km = haversine_km(float(r['lon']), float(r['lat']))
        print('  burn %d  (panels %s)' % (i, 'c,d' if i == 1 else 'e,f'))
        print('    event_id   %s' % r['event_id'])
        print('    date       %s' % r['burn_date'])
        print('    province   %s' % repair(r['province']))
        print('    area       %s ha' % r['area_ha'])
        print('    pre-NDVI   %.3f    slope %.2f deg' % (float(r['pre_ndvi']),
                                                         float(r['slope_mean'])))
        print('    distance   %.0f km from the rupture point' % km)

    d = np.array([haversine_km(float(r['lon']), float(r['lat'])) for r in rows])
    a = np.array([float(r['area_ha']) for r in rows])
    print('\nSECTION 4.5, distances from the rupture point')
    print('  median      %.0f km' % np.median(d))
    print('  IQR         %.0f - %.0f km' % (np.percentile(d, 25),
                                            np.percentile(d, 75)))
    print('  within 500  %d of %d (%.0f %%)' % ((d <= 500).sum(), len(d),
                                                100 * (d <= 500).mean()))
    print('  min, max    %.0f, %.0f km' % (d.min(), d.max()))
    print('  area_ha     median %.0f, IQR %.0f-%.0f'
          % (np.median(a), np.percentile(a, 25), np.percentile(a, 75)))


if __name__ == '__main__':
    main()
