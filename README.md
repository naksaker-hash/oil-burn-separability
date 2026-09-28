# oil-burn-separability

An open workflow for quantifying spectral confusion between land surface
classes, applied to terrestrial oil spills and fire scars in Sentinel-2
imagery.

This repository holds the code and the archived sample table behind the
manuscript submitted to *Computers & Geosciences*. Every number, figure
and table in that manuscript and its supplementary material can be
regenerated from what is here.

## What the code does

Remote sensing of a rare land surface target is limited less by
detection than by confusion with surfaces that look like the target.
This code measures that confusion. It takes a table of per-pixel
spectral features labelled by class and by event, and it reports how
separable each class pair is, with an uncertainty interval that reflects
the number of independent events rather than the number of pixels.

Four choices distinguish it from a conventional separability study.

- **Direction-free separability.** Separability is computed as
  `S = |AUC - 0.5| * 2` from the rank-sum identity with mid-ranks for
  ties. Folding at chance means a feature that responds in the
  unexpected direction is scored on its information content rather than
  on its sign.
- **Cluster bootstrap over events.** The bootstrap resamples whole
  events with replacement and redraws every pixel belonging to a
  selected event. Pixel-level resampling treats spatially dependent
  pixels as independent evidence and produces intervals that are far too
  narrow.
- **Nested leave-one-event-out selection.** Feature selection is
  repeated inside the training partition of every fold, so the reported
  generalisation figure carries no selection leakage and the stability
  of the selected pair becomes a measured quantity.
- **Label independence.** Class labels are delineated independently of
  every feature under evaluation, and label quality is then checked with
  a per-event statistic computed on a quantity that is measured but
  never selected upon.

The components are not specific to oil and fire and apply to any pair of
spectrally similar targets.

## Requirements

Python 3.9 or later, and two packages.

```
numpy>=1.24
matplotlib>=3.6
```

Install them with:

```
pip install -r requirements.txt
```

That is the whole dependency set. The rank-sum AUC, the
Jeffries-Matusita distance and the regularised Fisher linear
discriminant are implemented directly rather than taken from a
statistics library, so no compiled statistics or dataframe package is
needed. Re-running the sampling step additionally requires a Google
Earth Engine account, but the archived sample table is included here so
that the analysis reproduces without one.

No graphics processing unit is required. The complete analysis runs in
under ten minutes on a standard laptop with 8 GB of memory.

## Quick test

The repository includes a generator for a small synthetic sample with
the same schema as the real data, so the whole analysis path can be
exercised before touching anything real.

```
python3 make_test_sample.py test_sample.csv
python3 oil_burn_separability.py test_sample.csv
python3 nested_loeo.py test_sample.csv
python3 sequential_screen.py test_sample.csv
```

The first command writes `test_sample.csv`, which holds 150 synthetic
pixels for each of 11 events across 5 classes. The remaining commands
should each run to completion in a few seconds and write into
`results/` and `figures/`.

The synthetic values are drawn from normal distributions with a
deliberate offset between classes. **They are not real observations and
must not be used for any scientific claim.** Their only purpose is to
prove the code runs end to end.

## Reproducing the published results

Run the modules against the archived sample table in `data/`.

```
python3 oil_burn_separability.py data/oil_burn_samples.csv
python3 nested_loeo.py data/oil_burn_samples.csv
python3 sequential_screen.py data/oil_burn_samples.csv
python3 oil_mask_sensitivity.py
python3 label_quality.py
python3 temporal_lags.py
python3 sens_pre2018.py data/oil_burn_samples.csv
python3 sens_2024.py data/oil_burn_samples.csv
python3 make_supp_table.py
python3 make_supp_tables.py
```

Each module reads the same sample table and writes plain
comma-separated output, so intermediate results can be inspected without
rerunning upstream steps.

| Module | Produces |
| --- | --- |
| `oil_burn_separability.py` | The full 50-feature separability table, the dissociation figure, the per-event figure and the feature-space figure |
| `nested_loeo.py` | The nested leave-one-fire-out selection and validation |
| `sequential_screen.py` | The two-stage detect-then-reject screen at operational prevalence |
| `oil_mask_sensitivity.py` | Sensitivity of the result to the plume boundary |
| `label_quality.py` | The per-event label-quality statistic |
| `temporal_lags.py` | Acquisition lag summaries |
| `sens_pre2018.py`, `sens_2024.py` | Temporal subsetting checks |
| `make_supp_table.py`, `make_supp_tables.py` | Supplementary Tables S1 to S7 |
| `gee/narli_burn_confusion_gee.js` | The Earth Engine script that produced the sample table |

The random seeds used for sampling, for the cluster bootstrap and for
the cross-validation splits are fixed in source.

## A note on reproducibility

Earth Engine's stratified sampling is not bit-reproducible between
executions. Repeated runs of an identical script return slightly
different pixel counts. The archived extract at `data/oil_burn_samples.csv`, rather than a
re-execution of the Earth Engine script, is therefore the data of record
for every published value. The Python analysis run against that extract
is deterministic.

## Data sources

| Source | Identifier |
| --- | --- |
| Sentinel-2 L2A surface reflectance | `COPERNICUS/S2_SR_HARMONIZED` |
| Cloud and shadow screening | `GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED` |
| Copernicus GLO-30 DEM | `COPERNICUS/DEM/GLO30`, doi 10.5270/ESA-c5d3d65 |
| OpenLandMap soil texture class | `OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02`, doi 10.5281/zenodo.1475451 |
| EFFIS burnt area perimeters | forest-fire.emergency.copernicus.eu, downloaded 11 September 2026 |
| MODIS burned area, diagnostic only | `MODIS/061/MCD64A1` |

## Licence

MIT. See `LICENSE`.

## Citation

Aksaker, N. An open workflow for quantifying spectral confusion between
land surface classes, with application to terrestrial oil spills and
fire scars in Sentinel-2 imagery. Submitted to *Computers &
Geosciences*.

## Contact

Nazım Aksaker, Space Science and Solar Energy Research and Application
Center (UZAYMER), Çukurova University, Adana, Türkiye.
naksaker@cu.edu.tr
