#!/usr/bin/env python3
"""
Figure 1: where the spill is.

Two panels at two scales, in the visual idiom of Figure 5 so the two
maps read as a pair.

  (a) Türkiye, with the rupture point marked.
  (b) The 3 x 3 km analysis window in Sentinel-2 true colour on the
      first post-event acquisition, with the recovered plume outlined.

Panel (b) is rendered from oil_window_pixels.csv, the same export that
relocated the sampled pixels, so the image and the labels come from one
source and cannot disagree. Panel (a) needs a country boundary; run
fig01_outline_gee.js to produce it. Without it the panel is left blank
with a note rather than drawn from some other boundary, because Figure 5
uses GAUL and two maps in one paper should not use two coastlines.

Usage:
    python3 make_fig01_site.py
"""

import csv
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import utm

WINDOW_CSV = 'oil_window_pixels.csv'
PLUME_JSON = 'oil_plume_polygon.geojson'
OUTLINE = 'fig01_turkey_outline.geojson'

# Kirkuk--Ceyhan pipeline route, from OpenStreetMap via the Overpass
# API (c) OpenStreetMap contributors, ODbL. The segments were selected
# by connectivity from the one BOTAS oil pipeline way passing within
# 1.3 km of the rupture point, which identifies the line without
# relying on its name; the Baku--Tbilisi--Ceyhan line, which also
# terminates at Ceyhan, was excluded by name and lies 71 km away at its
# closest approach. OSM geometry is approximate and is drawn here to
# show the route, not to locate the pipe.
PIPELINE = 'fig01_pipeline_osm.geojson'

# World boundaries for the context panel, same GAUL source as the
# Türkiye outline. Run fig01_world_gee.js to produce it.
WORLD = 'fig01_world_outline.geojson'

# Antarctica and the far south are cropped: they carry no information
# here and squaring the panel off to them would halve the size of
# everything else.
WORLD_LAT = (-58.0, 84.0)
OUT = os.path.join('figures', 'fig01_site.png')

SPILL_LON, SPILL_LAT = 37.130123, 37.350183
ZONE = 37
PIXEL_M = 20.0

# Percentile stretch for the true-colour panel. Applied to all three
# bands together rather than per band, so the rendering does not
# colour-balance the plume away.
STRETCH_LO, STRETCH_HI = 2.0, 98.0

# Panel (b) crop, metres from the rupture point. The full window is
# 3 x 3 km, most of which is context the figure does not need; cropping
# makes the plume legible at column width and gives the panel an aspect
# closer to the locator beside it. The whole window is still shown in
# the pixel export and in Figure 4.
CROP_E, CROP_W = 1050.0, 700.0
CROP_N, CROP_S = 680.0, 680.0


def rings(geom):
    t = geom.get('type')
    if t == 'Polygon':
        for r in geom['coordinates']:
            yield np.asarray(r, dtype=float)
    elif t == 'MultiPolygon':
        for poly in geom['coordinates']:
            for r in poly:
                yield np.asarray(r, dtype=float)
    elif t == 'GeometryCollection':
        for g in geom.get('geometries', []):
            for r in rings(g):
                yield r


def load_window():
    with open(WINDOW_CSV, newline='') as fh:
        rows = list(csv.DictReader(fh))
    xy = [utm.forward(float(r['longitude']), float(r['latitude']), ZONE)
          for r in rows]
    xs = np.array([p[0] for p in xy])
    ys = np.array([p[1] for p in xy])
    x0, y0 = xs.min(), ys.min()
    i = np.round((xs - x0) / PIXEL_M).astype(int)
    j = np.round((ys - y0) / PIXEL_M).astype(int)
    img = np.full((j.max() + 1, i.max() + 1, 3), np.nan)
    for k, r in enumerate(rows):
        img[j[k], i[k]] = [float(r['B4']), float(r['B3']), float(r['B2'])]
    return img, x0, y0


def main():
    os.makedirs('figures', exist_ok=True)
    if not os.path.exists(WINDOW_CSV):
        raise SystemExit('Not found: %s' % WINDOW_CSV)

    img, x0, y0 = load_window()
    lo = np.nanpercentile(img, STRETCH_LO)
    hi = np.nanpercentile(img, STRETCH_HI)
    disp = np.nan_to_num(np.clip((img - lo) / (hi - lo), 0, 1), nan=1.0)
    H, W = img.shape[0], img.shape[1]

    ox, oy = utm.forward(SPILL_LON, SPILL_LAT, ZONE)
    ex = [-(ox - x0), W * PIXEL_M - (ox - x0),
          -(oy - y0), H * PIXEL_M - (oy - y0)]

    # Panels run coarse to fine, left to right: world, country, site.
    fig, (ax_world, ax_country, ax_site) = plt.subplots(
        1, 3, figsize=(13, 3.3),
        gridspec_kw={'width_ratios': [1.15, 1.35, 1.0]})

    # ---- panel (a) locator ----
    if os.path.exists(OUTLINE):
        with open(OUTLINE) as fh:
            gj = json.load(fh)
        for ft in gj.get('features', [gj]):
            for r in rings(ft.get('geometry', ft)):
                ax_country.fill(r[:, 0], r[:, 1], facecolor='#f2f2f2',
                         edgecolor='#999999', lw=0.7, zorder=1)
        ax_country.set_aspect(1.0 / np.cos(np.deg2rad(SPILL_LAT)))

        # Country bounds fix the frame, so the pipeline running on into
        # Iraq is clipped at the edge rather than stretching the panel.
        xs = np.concatenate([r[:, 0] for ft in gj.get('features', [gj])
                             for r in rings(ft.get('geometry', ft))])
        ys = np.concatenate([r[:, 1] for ft in gj.get('features', [gj])
                             for r in rings(ft.get('geometry', ft))])
        ax_country.set_xlim(xs.min() - 0.4, xs.max() + 0.4)
        ax_country.set_ylim(ys.min() - 0.3, ys.max() + 0.3)
    else:
        ax_country.text(0.5, 0.5, 'run fig01_outline_gee.js\nfor the boundary',
                 ha='center', va='center', fontsize=9, color='#999999',
                 transform=ax_country.transAxes)
        print('NOTE: %s missing; panel (a) left blank.' % OUTLINE)

    if os.path.exists(PIPELINE):
        with open(PIPELINE) as fh:
            pl = json.load(fh)
        for k, ft in enumerate(pl['features']):
            c = np.asarray(ft['geometry']['coordinates'], dtype=float)
            ax_country.plot(c[:, 0], c[:, 1], c='#c0392b', lw=1.3, zorder=2,
                     solid_capstyle='round', clip_on=True,
                     label='Kirkuk\u2013Ceyhan pipeline' if k == 0 else None)
        # Both ends are labelled "to", not named: OpenStreetMap's
        # coverage stops short of the Ceyhan terminal, and the Kirkuk
        # end lies outside the frame. Naming the drawn endpoints would
        # claim more than the data supports.
        pts = np.vstack([np.asarray(ft['geometry']['coordinates'],
                                    dtype=float)
                         for ft in pl['features']])
        ylo, yhi = ax_country.get_ylim()
        west = pts[np.argmin(pts[:, 0])]
        vis = pts[(pts[:, 1] >= ylo) & (pts[:, 1] <= yhi)]
        south = vis[np.argmin(vis[:, 1])] if len(vis) else None
        ax_country.annotate('to Ceyhan', west, xytext=(-5, 9),
                     textcoords='offset points', fontsize=8, ha='right')
        if south is not None:
            ax_country.annotate('to Kirkuk', south, xytext=(7, 4),
                         textcoords='offset points', fontsize=8)
    else:
        print('NOTE: %s missing; pipeline not drawn.' % PIPELINE)

    ax_country.scatter([SPILL_LON], [SPILL_LAT], s=200, marker='*',
                facecolor='#f1c40f', edgecolor='black', lw=0.8, zorder=4)
    ax_country.annotate('Narlı', (SPILL_LON, SPILL_LAT),
                 xytext=(8, -12), textcoords='offset points', fontsize=9)
    ax_country.set_xlabel('longitude (°E)', fontsize=9)
    ax_country.set_ylabel('latitude (°N)', fontsize=9)
    ax_country.set_title('(b) Kirkuk–Ceyhan pipeline rupture, 19 June 2018',
                  fontsize=9, loc='left')
    if os.path.exists(PIPELINE):
        ax_country.legend(loc='upper left', frameon=False, fontsize=8)
    ax_country.tick_params(labelsize=8)
    ax_country.grid(True, lw=0.3, c='#dddddd', zorder=0)

    # ---- panel (b) the site ----
    ax_site.imshow(np.flipud(disp), origin='upper', extent=ex,
               interpolation='nearest')
    ax_site.set_xlim(-CROP_W, CROP_E)
    ax_site.set_ylim(-CROP_S, CROP_N)

    if os.path.exists(PLUME_JSON):
        with open(PLUME_JSON) as fh:
            pg = json.load(fh)
        for r in rings(pg['geometry']):
            pxy = np.array([utm.forward(a, b, ZONE) for a, b in r])
            ax_site.plot(pxy[:, 0] - ox, pxy[:, 1] - oy, c='#c0392b', lw=1.4,
                     zorder=3)
        lab = 'digitised plume, 8.72 ha'
    else:
        lab = None
        print('NOTE: %s missing; plume outline omitted.' % PLUME_JSON)

    ax_site.scatter([0], [0], s=200, marker='*', facecolor='#f1c40f',
                edgecolor='black', lw=0.8, zorder=4)

    # Scale bar, sized from the panel rather than fixed, so it stays
    # correct if the window changes.
    bar = 250.0
    bx0 = -CROP_W + 0.07 * (CROP_E + CROP_W)
    by = -CROP_S + 0.09 * (CROP_N + CROP_S)
    ax_site.plot([bx0, bx0 + bar], [by, by], c='white', lw=3.2, zorder=5)
    ax_site.plot([bx0, bx0 + bar], [by, by], c='black', lw=1.6, zorder=6)
    ax_site.text(bx0 + bar / 2, by + 0.03 * (CROP_N + CROP_S), '250 m',
             ha='center', va='bottom', fontsize=8, color='black', zorder=6)

    handles = [Line2D([], [], ls='', marker='*', mfc='#f1c40f',
                      mec='black', ms=13, label='rupture point')]
    if lab:
        handles.append(Line2D([], [], c='#c0392b', lw=1.4, label=lab))
    ax_site.legend(handles=handles, loc='upper left', frameon=True,
               framealpha=0.85, fontsize=8)

    ax_site.set_xlabel('metres east of the rupture point', fontsize=9)
    ax_site.set_ylabel('metres north', fontsize=9)
    ax_site.set_aspect('equal')
    ax_site.set_title('(c) Sentinel-2 true colour, 20 June 2018', fontsize=9,
                  loc='left')
    ax_site.tick_params(labelsize=8)

    # ---- panel (c) world context ----
    if os.path.exists(WORLD):
        with open(WORLD) as fh:
            wj = json.load(fh)
        for ft in wj.get('features', [wj]):
            name = (ft.get('properties') or {}).get('ADM0_NAME', '')
            hit = (name == 'Turkey')
            for r in rings(ft.get('geometry', ft)):
                ax_world.fill(r[:, 0], r[:, 1],
                         facecolor='#c0392b' if hit else '#ececec',
                         edgecolor='#c0392b' if hit else '#bbbbbb',
                         lw=0.6 if hit else 0.25,
                         zorder=3 if hit else 1)
        ax_world.set_xlim(-180, 180)
        ax_world.set_ylim(*WORLD_LAT)
        ax_world.set_aspect(1.0 / np.cos(np.deg2rad(39.0)))
        ax_world.annotate('T\u00fcrkiye', (SPILL_LON, SPILL_LAT),
                     xytext=(10, 26), textcoords='offset points',
                     fontsize=8, color='#c0392b',
                     arrowprops=dict(arrowstyle='-', lw=0.7,
                                     color='#c0392b'))
    else:
        ax_world.text(0.5, 0.5, 'run fig01_world_gee.js\nfor the world outline',
                 ha='center', va='center', fontsize=9, color='#999999',
                 transform=ax_world.transAxes)
        print('NOTE: %s missing; panel (c) left blank.' % WORLD)

    ax_world.set_title('(a) Regional context', fontsize=9, loc='left')
    ax_world.set_xticks([])
    ax_world.set_yticks([])
    for sp in ax_world.spines.values():
        sp.set_linewidth(0.6)

    fig.tight_layout()
    fig.savefig(OUT, dpi=300)
    plt.close(fig)
    print('written: %s  (%d x %d px panel)' % (OUT, W, H))


if __name__ == '__main__':
    main()
