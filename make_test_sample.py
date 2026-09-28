#!/usr/bin/env python3
"""
Generate a small synthetic sample table with the same schema as
oil_burn_samples.csv.

This exists so that the analysis path can be exercised without the real
data, which is the quick test required by the Computers & Geosciences
repository guidelines. The numbers it produces are drawn from simple
normal distributions with a deliberate offset between classes. They are
not real observations and must not be used for any scientific claim.
They exist only to prove that the code runs end to end.

Usage:
    python3 make_test_sample.py [output.csv] [n_per_event]

Default output is test_sample.csv with 150 pixels per event. The
analysis modules require at least MIN_EVENT_N = 100 pixels before an
event is admitted to the event-level bootstrap and the
leave-one-event-out loop, so smaller values silently produce empty
result tables.
"""

import csv
import sys

import numpy as np

# Schema of oil_burn_samples.csv, in order.
META = ['class_name', 'event_id', 'event_date', 'event_end_date',
        'pre_date', 'post_date', 'area_ha', 'pre_ndvi', 'agri_pct',
        'slope_mean', 'soil_mode', 'province', 'pre_NDVI_px',
        'slope_px', 'soil_px']

BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']

INDICES = ['NBR', 'NBR2', 'BAI', 'logBAI', 'MIRBI', 'BAIS2', 'FMR',
           'NDBI', 'ALBEDO', 'NDVI', 'DVI', 'ARVI', 'ARI', 'IRECI',
           'NDWI']

# The differenced features carry a d_ prefix and cover both the bands
# and the indices.
DIFFS = ['d_' + f for f in BANDS + INDICES]

COLUMNS = META + BANDS + INDICES + DIFFS

# Class offsets. Oil and burn are placed close together and clean soil
# further away, which reproduces the qualitative structure the paper
# reports without reproducing any actual value.
OFFSET = {'oil': 0.00, 'oil_edge': 0.05, 'burn': 0.15,
          'control_clean': 1.20, 'control_unburned': 1.40}

# Synthetic events. The burn side needs several events for the cluster
# bootstrap and the leave-one-event-out loop to have anything to
# resample and withhold. Class names match those the analysis modules
# select on.
EVENTS = ([('oil', 'OIL_TEST_0001'), ('oil_edge', 'OIL_TEST_0001')]
          + [('burn', 'BURN_TEST_%04d' % i) for i in range(1, 9)]
          + [('control_clean', 'CLEAN_TEST_0001'),
             ('control_unburned', 'CONTROL_TEST_0001')])


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else 'test_sample.csv'
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 150

    rng = np.random.default_rng(20260928)

    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(COLUMNS)

        for cls, event_id in EVENTS:
            base = OFFSET[cls]
            # Give each event its own small shift so that the
            # between-event variance the cluster bootstrap is meant to
            # capture actually exists in the test data.
            shift = rng.normal(0.0, 0.05)

            for _ in range(n):
                meta = [
                    cls,
                    event_id,
                    '2018-06-15',
                    '2018-06-15',
                    '2018-06-10',
                    '2018-06-20',
                    '%.2f' % rng.uniform(1.0, 50.0),
                    '%.3f' % rng.uniform(0.10, 0.35),
                    '%.1f' % rng.uniform(0.0, 100.0),
                    '%.2f' % rng.uniform(0.5, 8.0),
                    str(int(rng.integers(1, 12))),
                    'TEST',
                    '%.3f' % rng.uniform(0.10, 0.35),
                    '%.2f' % rng.uniform(0.5, 8.0),
                    str(int(rng.integers(1, 12))),
                ]
                vals = rng.normal(base + shift, 0.30,
                                  len(BANDS) + len(INDICES) + len(DIFFS))
                w.writerow(meta + ['%.6f' % v for v in vals])

    rows = n * len(EVENTS)
    print('wrote %s' % out)
    print('  %d rows, %d columns' % (rows, len(COLUMNS)))
    print('  %d events across %d classes'
          % (len(EVENTS), len(set(c for c, _ in EVENTS))))
    print('')
    print('This is synthetic data. Do not use it for any scientific claim.')


if __name__ == '__main__':
    main()
