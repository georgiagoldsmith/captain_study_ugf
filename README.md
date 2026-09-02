# Conservation prioritisation in the Upper Guinean Forest

Compares two approaches to choosing where to protect land: **CAPTAIN**, a
reinforcement-learning model that simulates an ecosystem responding to
protection over time, and **prioritizr**, which solves the same problem exactly
as an integer program against fixed inputs.

Both select roughly 30% of the study area, cover the same 18 bird species, and
run on the same 3 km grid, so their solutions are directly comparable.

## What is and is not in this repository

The code is here. **Most of the data is not** — `data/` is around 130 MB of
rasters, and `data/gHM/gHM.tif` alone is 2.1 GB. Those are ignored by git.

Tracked: the study boundary, the species trait tables, the protected-area
clip, and `data/captain/disturbance.tif` (the grid template). These are small
and awkward to rebuild.

Not tracked, and needed before anything will run:

| Path | What | Source |
|---|---|---|
| `data/gHM/gHM.tif` | present-day human modification | Theobald et al. 2025 |
| `data/birds_aoh/clipped/` | 18 species Area of Habitat rasters | Lumbierres et al. 2022 (Dryad) |
| `data/cocoa suitability/` | cocoa crop suitability | CropSuite v1.0, Zabel et al. 2024 |
| `data/future lulc/SSP4_RCP34/` | projected land cover, 2020 and 2075 | Zhang et al. 2023 |
| `data/LC/` | land cover clips (urban, cropland, LCCS) | C3S LCCS 300 m |

`outputs/` is ignored entirely: everything in it is produced by the scripts.

## Running it

CAPTAIN runs need the `captain3preview` Python environment and a separate data
directory outside this repo; prioritizr needs Gurobi.

```
# 1. prepare inputs        (only if the derived layers need rebuilding)
Rscript scripts/spatial_data_prep.R

# 2. prioritizr            writes outputs/prioritizr/ and outputs/layers_3km/
Rscript scripts/prioritizr_ugf_prioritization.R

# 3. CAPTAIN               export inputs, then train and predict
Rscript scripts/captain_data_prep.R
Rscript scripts/create_future_disturbance.R
python  scripts/captain_train.py
python  scripts/captain_predict.py
bash    scripts/run_future_staged.sh     # the staged-protection variant

# 4. results
Rscript scripts/captain_seed_comparison.R
Rscript scripts/plot_captain_figures.R
Rscript scripts/agreement_maps_captain_prioritizr.R
Rscript scripts/compare_characteristics_captain_prioritizr.R
```

Step 2 must precede step 3: `captain_data_prep.R` reads the layers prioritizr
writes, so both models are fed identical numbers.

## Scripts

| Script | Does |
|---|---|
| `spatial_data_prep.R` | clips land cover, cocoa and protected areas to the UGF |
| `ugf_birds_iucn_and_aoh.R` | IUCN status per species, matches them to AOH rasters |
| `carrying_capacity_estimation.R` | per-species carrying capacity from density and body mass |
| `prioritizr_ugf_prioritization.R` | the prioritizr solution, its coverage table and map |
| `captain_data_prep.R` | writes CAPTAIN's input layers on the 3 km grid |
| `create_future_disturbance.R` | 2075 disturbance from projected land cover |
| `captain_train.py` / `captain_predict.py` | train a policy, then turn it into a plan |
| `run_future_staged.sh` | the staged-protection run (161 cells per step) |
| `captain_seed_comparison.R` | compares the 10 seeds |
| `plot_captain_figures.R` | the CAPTAIN solution maps |
| `agreement_maps_captain_prioritizr.R` | where the two models diverge |
| `compare_characteristics_captain_prioritizr.R` | what each model selected |

## Known gaps

`carrying_capacity_estimation.R` reads TetraDENSITY and AVONET from `/tmp`,
which no longer exist. Its outputs are committed, so nothing downstream breaks,
but the step cannot currently be re-run without re-downloading both.

`spatial_data_prep.R` and `ugf_birds_iucn_and_aoh.R` also depend on source data
not kept here (WDPA source shapefiles, OSM roads, the raw AOH order folders).
They are provenance records rather than steps you re-run.
