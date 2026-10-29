# =============================================================================
# Nagaon watershed cLHS study: run the steps one at a time in R
# Open R (R.app), then copy each line below into the console.
# =============================================================================

# Set the working directory to the folder that contains this file, e.g.
# setwd("~/Documents/nagaon-clhs")

source("R/00_setup.R")        # once: installs packages and WhiteboxTools (5-10 min)

source("R/01_watershed.R")    # downloads boundaries and DEM, lists candidate watersheds (10-20 min)
#   -> open outputs/figures/01_watershed_candidates.png and outputs/watershed_candidates.csv
#   -> choose a watershed and put its rank number in WS_RANK at the top of R/02_covariates.R

source("R/02_covariates.R")   # terrain, Sentinel-2, land cover, flooding, roads (15-30 min)
source("R/03_clhs.R")         # selects 30 profile sites + backups (5-10 min)
source("R/04_figures.R")      # maps and charts for the presentation (1 min)

# Results:
#   outputs/field/nagaon_ws_clhs_sites.csv / .kml / .gpx   field points
#   outputs/figures/                                         figures for the slides
#   outputs/clhs_diagnostics.csv, outputs/clhs_vs_random.csv tables
