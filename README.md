# Open-data R workflow for selecting soil profile sites by conditioned Latin hypercube sampling (cLHS)

R scripts accompanying:

> Tamuly, D. and Dutta, S. An open-data R workflow for selecting soil profile sites by conditioned Latin hypercube sampling in a hill–floodplain watershed of Nagaon, Assam. *Journal of Soil and Water Conservation* (submitted).

The workflow selects a sub-watershed, builds a 30 m covariate stack from public data, applies cLHS with access, spacing and replacement-site constraints, and writes field-ready site lists (CSV, KML, GPX, GeoPackage). Everything runs in R; no licensed data or desktop GIS steps are needed. It can be applied to any district by changing the district name in `01_watershed.R`.

## Steps

| Script | What it does | Main outputs |
|---|---|---|
| `00_setup.R` | Installs R packages and WhiteboxTools; creates folders | – |
| `01_watershed.R` | District boundary (GADM 4.1), HydroBASINS level 9–10 sub-watersheds, Copernicus GLO-30 DEM; ranks candidates by relief | `outputs/watershed_candidates.csv`, figure |
| `02_covariates.R` | Terrain (elevation, slope, TPI), hydrology (TWI, HAND, distance to stream), Sentinel-2 dry-season indices (NDVI, BSI, MNDWI), ESA WorldCover, JRC surface water, OpenStreetMap roads; non-soil mask | `data/*.tif` |
| `03_clhs.R` | Candidate points, 1 km road-access filter, cLHS (`clhs` package, 5 runs × 50,000 iterations), 1.5 km minimum spacing, up to three replacement sites per profile, Monte Carlo comparison with 500 random designs | `outputs/field/*`, `outputs/clhs_diagnostics.csv`, `outputs/clhs_vs_random.csv` |
| `04_figures.R` | Maps and diagnostic plots | `outputs/figures/*.png` |

Run them in order from `RUN_ME.R`, with the working directory set to the folder that holds the scripts. Set `WS_RANK` at the top of `02_covariates.R` to the sub-watershed chosen in step 1 (rank 1 in the paper). Key settings (`ACCESS_M`, `MIN_KM`, `ITER`, number of sites) are at the top of `03_clhs.R`. Random seeds are fixed (`set.seed(2026)`).

Downloaded data are cached in `data/raw/`, so later runs are faster. An internet connection is needed for the first run.

## Results from the paper

The three CSV files hold the outputs reported in the paper:

- `nagaon_ws_clhs_sites.csv` – 30 profile sites (NGW-P01 to P30) and 86 replacement sites, with WGS 84 coordinates and covariate values
- `clhs_diagnostics.csv` – stratum coverage and Kolmogorov–Smirnov statistics per covariate
- `clhs_vs_random.csv` – cLHS compared with 500 simple random designs

## Requirements

R ≥ 4.3 with terra, sf, dplyr, tidyr, ggplot2, patchwork, clhs (0.9.2 used), whitebox, rstac, osmdata and geodata. Tested with R 4.5.2 on macOS.

## Data sources

Copernicus GLO-30 DEM (ESA); Sentinel-2 L2A (ESA/Copernicus, via Earth Search STAC); ESA WorldCover 2021; JRC Global Surface Water; HydroBASINS (HydroSHEDS); GADM 4.1; OpenStreetMap contributors (ODbL).

## Citation

Please cite the paper above and this archive (see `CITATION.cff`).

## Licence

Code: MIT licence (see `LICENSE`). Site list and diagnostics (CSV files): CC BY 4.0.

## Contact

Dr. Danish Tamuly, Department of Soil Science, Assam Agricultural University, Jorhat 785013, Assam, India – danish.tamuly@aau.ac.in
