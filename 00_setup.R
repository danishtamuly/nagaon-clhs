# =============================================================================
# Nagaon watershed cLHS study — Step 0: set up R
# Dr. Danish Tamuly, Department of Soil Science, AAU, Jorhat
# Run once. Installs the packages used in steps 1-4.
# =============================================================================

pkgs <- c(
  "terra",      # raster processing
  "sf",         # vector data (boundaries, points)
  "dplyr",      # tables
  "ggplot2",    # figures
  "tidyr",
  "patchwork",  # combine figures
  "clhs",       # conditioned Latin hypercube sampling (Roudier)
  "whitebox",   # hydrology: flow, streams, HAND, TWI
  "rstac",      # find Sentinel-2 scenes
  "osmdata",    # roads from OpenStreetMap
  "geodata"     # administrative boundaries (GADM)
)
new <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

# WhiteboxTools engine (downloaded once, ~30 MB)
library(whitebox)
if (!wbt_init()) whitebox::install_whitebox()

# Project folders
for (d in c("data", "data/raw", "outputs", "outputs/figures", "outputs/field"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("Setup complete. R", R.version.string, "\n")
