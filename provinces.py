#!/usr/bin/env python3
"""
Repair Turkish province names in oil_burn_samples.csv.

The `province` column is not merely mis-decoded, it is damaged. The
non-ASCII letters were replaced before the file was written: most became
a literal question mark (U+003F), and two became the Unicode replacement
character (U+FFFD). Neither is recoverable by reading the file in a
different encoding, because the original bytes are gone.

They are recoverable by name, because Turkey has 81 provinces and each
damaged string matches exactly one of them. That is what this module
does. It is a lookup, not a decode, and it is deliberately strict: an
unrecognised string raises rather than passing through, so a future
export that breaks in a new way fails loudly instead of quietly
producing a 24th province.

The column arrives from the EFFIS inventory asset via
`feat.get('PROVINCE')` in narli_burn_confusion_gee.js, which only passes
it through. If the export is ever re-run, normalise at that source.

Usage:
    python3 provinces.py [path/to/oil_burn_samples.csv]

    from provinces import repair
    name = repair(row['province'])
"""

import collections
import csv
import os
import sys

# Damaged form -> correct name. Every key was observed in the archived
# oil_burn_samples.csv; the values are the standard Turkish spellings.
#
# Two damage patterns appear, which is why some keys carry '?' and
# others U+FFFD: the letters absent from Latin-1 (S-cedilla, dotless i,
# g-breve) became '?', while two letters present in Latin-1 (o-diaeresis,
# C-cedilla) became U+FFFD. The cause is upstream and unknown; the
# mapping is unaffected by it.
REPAIR = {
    '?anl?urfa':   'Şanlıurfa',
    'Diyarbak?r':  'Diyarbakır',
    'Ad?yaman':    'Adıyaman',
    'Elaz??':      'Elazığ',
    'Mu?':         'Muş',
    'Tekirda?':    'Tekirdağ',
    '??rnak':      'Şırnak',
    '�orum':  'Çorum',
    'Bing�l': 'Bingöl',
}

# Names that survived intact, listed so that `repair` can tell a clean
# value from an unrecognised one rather than assuming anything ASCII is
# fine.
INTACT = {
    'Adana', 'Aksaray', 'Ankara', 'Antalya', 'Batman', 'Gaziantep',
    'Kayseri', 'Kilis', 'Konya', 'Malatya', 'Mardin', 'Siirt', 'Sivas',
    'Yozgat',
}


def repair(name):
    """Correct one province string. Raises on anything unrecognised."""
    if name in INTACT:
        return name
    if name in REPAIR:
        return REPAIR[name]
    raise ValueError(
        'unrecognised province %r: the export has damaged a name in a '
        'way this map does not cover. Add it to REPAIR rather than '
        'letting it through.' % name)


def repaired_column(rows, field='province'):
    return [repair(r[field]) for r in rows]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'oil_burn_samples.csv'
    if not os.path.exists(path):
        raise SystemExit('Not found: %s' % path)
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))

    burn = [r for r in rows if r['class_name'] == 'burn']

    # Provinces are a property of the event, not the pixel, so count
    # events. Counting pixels would weight by fire size.
    ev = {}
    for r in burn:
        ev[r['event_id']] = repair(r['province'])
    counts = collections.Counter(ev.values())

    print('MATCHED BURN EVENTS BY PROVINCE  (n = %d events)' % len(ev))
    print('-' * 40)
    for name, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        print('%-14s %3d' % (name, n))
    print('-' * 40)
    print('%-14s %3d' % ('provinces', len(counts)))

    top4 = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:4]
    share = sum(n for _, n in top4) / float(len(ev))
    print()
    print('top four %s = %d of %d events (%.0f %%)'
          % (', '.join('%s %d' % (k, v) for k, v in top4),
             sum(n for _, n in top4), len(ev), 100 * share))

    print()
    print('LaTeX-safe forms for the manuscript (UTF-8 source, or T1 + '
          'inputenc):')
    for name, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:6]:
        print('  %s (%d)' % (name, n))


if __name__ == '__main__':
    main()
