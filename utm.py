#!/usr/bin/env python3
"""
WGS84 <-> UTM, Snyder's series. Pure Python, no dependencies.

Needed because the Sentinel-2 sampling grid is UTM, not geographic. At
this latitude the UTM grid is rotated about 1.2 degrees from true north
by meridian convergence, so in longitude and latitude no two pixels in
the 3 km window share a coordinate and the grid cannot be reconstructed
there. Projected to UTM it is exactly regular at 20 m.

Accurate to well under a millimetre over a UTM zone.
"""

import math

A = 6378137.0                    # WGS84 semi-major axis
F = 1.0 / 298.257223563          # flattening
E2 = F * (2 - F)                 # first eccentricity squared
EP2 = E2 / (1 - E2)              # second eccentricity squared
K0 = 0.9996                      # UTM scale factor
E0 = 500000.0                    # false easting


def zone_from_lon(lon):
    return int((lon + 180) / 6) + 1


def _lon0(zone):
    return math.radians(zone * 6 - 183)


def forward(lon, lat, zone):
    """(lon, lat) in degrees -> (easting, northing) in metres."""
    p = math.radians(lat)
    l = math.radians(lon)
    l0 = _lon0(zone)

    N = A / math.sqrt(1 - E2 * math.sin(p) ** 2)
    T = math.tan(p) ** 2
    C = EP2 * math.cos(p) ** 2
    Aa = (l - l0) * math.cos(p)

    M = A * ((1 - E2 / 4 - 3 * E2 ** 2 / 64 - 5 * E2 ** 3 / 256) * p
             - (3 * E2 / 8 + 3 * E2 ** 2 / 32 + 45 * E2 ** 3 / 1024) * math.sin(2 * p)
             + (15 * E2 ** 2 / 256 + 45 * E2 ** 3 / 1024) * math.sin(4 * p)
             - (35 * E2 ** 3 / 3072) * math.sin(6 * p))

    x = K0 * N * (Aa + (1 - T + C) * Aa ** 3 / 6
                  + (5 - 18 * T + T ** 2 + 72 * C - 58 * EP2) * Aa ** 5 / 120) + E0
    y = K0 * (M + N * math.tan(p)
              * (Aa ** 2 / 2 + (5 - T + 9 * C + 4 * C ** 2) * Aa ** 4 / 24
                 + (61 - 58 * T + T ** 2 + 600 * C - 330 * EP2) * Aa ** 6 / 720))
    return x, y


def inverse(x, y, zone):
    """(easting, northing) in metres -> (lon, lat) in degrees."""
    l0 = _lon0(zone)
    M = y / K0
    mu = M / (A * (1 - E2 / 4 - 3 * E2 ** 2 / 64 - 5 * E2 ** 3 / 256))
    e1 = (1 - math.sqrt(1 - E2)) / (1 + math.sqrt(1 - E2))

    p1 = (mu
          + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * math.sin(2 * mu)
          + (21 * e1 ** 2 / 16 - 55 * e1 ** 4 / 32) * math.sin(4 * mu)
          + (151 * e1 ** 3 / 96) * math.sin(6 * mu)
          + (1097 * e1 ** 4 / 512) * math.sin(8 * mu))

    C1 = EP2 * math.cos(p1) ** 2
    T1 = math.tan(p1) ** 2
    N1 = A / math.sqrt(1 - E2 * math.sin(p1) ** 2)
    R1 = A * (1 - E2) / (1 - E2 * math.sin(p1) ** 2) ** 1.5
    D = (x - E0) / (N1 * K0)

    lat = p1 - (N1 * math.tan(p1) / R1) * (
        D ** 2 / 2
        - (5 + 3 * T1 + 10 * C1 - 4 * C1 ** 2 - 9 * EP2) * D ** 4 / 24
        + (61 + 90 * T1 + 298 * C1 + 45 * T1 ** 2 - 252 * EP2 - 3 * C1 ** 2) * D ** 6 / 720)
    lon = l0 + (D - (1 + 2 * T1 + C1) * D ** 3 / 6
                + (5 - 2 * C1 + 28 * T1 - 3 * C1 ** 2 + 8 * EP2 + 24 * T1 ** 2)
                * D ** 5 / 120) / math.cos(p1)
    return math.degrees(lon), math.degrees(lat)
