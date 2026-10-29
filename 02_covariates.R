# =============================================================================
# Step 2: Build the covariate stack (scorpan factors) for the chosen watershed
# Grid: 30 m, UTM zone 46N (EPSG:32646)
# r (relief)   : elevation, slope, TPI, TWI, HAND, distance to stream  [Copernicus GLO-30 + WhiteboxTools]
# o (organisms): dry-season NDVI, bare soil index, MNDWI                [Sentinel-2 L2A, Nov-Feb]
# o (land use) : ESA WorldCover 2021
# flooding     : JRC Global Surface Water occurrence 1984-2021
# access       : distance to road                                       [OpenStreetMap]
# =============================================================================
library(sf); library(terra); library(dplyr); library(ggplot2); library(whitebox)
sf_use_s2(FALSE)
options(timeout = 1200)
wbt_init()

# >>> Set the watershed you chose in step 1 (its rank number on the map) <<<
WS_RANK <- 1
CRS     <- "EPSG:32646"
RES     <- 30
wp <- function(f) file.path(getwd(), f)            # full paths for WhiteboxTools

ws <- st_read("data/watershed_candidates.gpkg", quiet = TRUE) |> filter(rank == WS_RANK)
stopifnot(nrow(ws) == 1)
WS_ID <- ws$HYBAS_ID
ws <- ws |> st_transform(32646)
st_write(ws, "data/watershed.gpkg", delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("Watershed %s: %.0f ha\n", WS_ID, as.numeric(st_area(ws)) / 1e4))

aoi  <- st_buffer(ws, 1500)                          # 1.5 km buffer for edge effects
aoiG <- st_transform(aoi, 4326)
e0   <- as.vector(ext(vect(aoi)))                   # snap the grid to whole 30 m cells
tmpl <- rast(xmin = floor(e0[1] / RES) * RES, xmax = ceiling(e0[2] / RES) * RES,
             ymin = floor(e0[3] / RES) * RES, ymax = ceiling(e0[4] / RES) * RES,
             resolution = RES, crs = CRS)
align <- function(r, method = "bilinear")           # put any layer exactly on the 30 m grid
  if (compareGeom(r, tmpl, stopOnError = FALSE)) r else resample(r, tmpl, method = method)

# 2a. DEM and terrain ------------------------------------------------------
dem <- rast("data/raw/dem_nagaon_cop30.tif") |> crop(vect(aoiG)) |> project(tmpl, method = "bilinear")
names(dem) <- "elev"
writeRaster(dem, "data/dem.tif", overwrite = TRUE)

slope <- terrain(dem, "slope", unit = "degrees"); names(slope) <- "slope"
w <- focalMat(dem, 300, "circle")                    # 300 m neighbourhood
tpi <- dem - focal(dem, w, fun = "sum", na.rm = TRUE); names(tpi) <- "tpi"

wbt_breach_depressions_least_cost(wp("data/dem.tif"), wp("data/dem_breached.tif"), dist = 50, fill = TRUE)
wbt_d8_flow_accumulation(wp("data/dem_breached.tif"), wp("data/fa_cells.tif"), out_type = "cells")
wbt_d8_flow_accumulation(wp("data/dem_breached.tif"), wp("data/sca.tif"), out_type = "specific contributing area")
wbt_slope(wp("data/dem_breached.tif"), wp("data/slope_wbt.tif"))
wbt_wetness_index(wp("data/sca.tif"), wp("data/slope_wbt.tif"), wp("data/twi.tif"))
stream_cells <- ceiling(1e6 / RES^2)                 # streams drain >= 1 km2
wbt_extract_streams(wp("data/fa_cells.tif"), wp("data/streams.tif"), threshold = stream_cells, zero_background = TRUE)
wbt_elevation_above_stream(wp("data/dem_breached.tif"), wp("data/streams.tif"), wp("data/hand.tif"))

twi  <- rast("data/twi.tif");  names(twi)  <- "twi"
hand <- rast("data/hand.tif"); names(hand) <- "hand"
streams <- rast("data/streams.tif")
dstream <- distance(classify(streams, cbind(0, NA))); names(dstream) <- "dist_stream"
hand <- classify(hand, cbind(-Inf, 0, 0))            # tiny negative values from DEM noise -> 0

# Cloud-hosted rasters are read with GDAL warp straight onto the 30 m grid.
# GDAL then uses the files' built-in overviews, which downloads far less data.
gdal_cfg <- c(GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
              CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.TIF",
              GDAL_HTTP_MAX_RETRY = "5", GDAL_HTTP_RETRY_DELAY = "5",
              GDAL_HTTP_MULTIRANGE = "YES", GDAL_HTTP_MERGE_CONSECUTIVE_RANGES = "YES",
              VSI_CACHE = "TRUE")
te <- as.vector(ext(tmpl))[c(1, 3, 2, 4)]
warp_to_grid <- function(url, out, method = "average", nodata = NULL) {
  if (file.exists(out)) return(rast(out))
  opt <- c("-t_srs", CRS, "-te", te, "-tr", RES, RES, "-r", method, "-ovr", "AUTO", "-overwrite")
  if (!is.null(nodata)) opt <- c(opt, "-srcnodata", nodata, "-dstnodata", nodata)
  tmp <- paste0(out, ".part.tif")
  sf::gdal_utils("warp", paste0("/vsicurl/", url), tmp, options = as.character(opt),
                 config_options = gdal_cfg, quiet = TRUE)
  file.rename(tmp, out)
  rast(out)
}
dir.create("data/raw/s2_30m", showWarnings = FALSE, recursive = TRUE)

# 2b. Sentinel-2 dry-season composite -------------------------------------
if (!file.exists("data/s2_composite.tif")) {
  library(rstac)
  items <- stac("https://earth-search.aws.element84.com/v1") |>
    stac_search(collections = "sentinel-2-c1-l2a", bbox = as.numeric(st_bbox(aoiG)),
                datetime = "2025-11-15T00:00:00Z/2026-02-28T23:59:59Z", limit = 100) |>
    ext_query(`eo:cloud_cover` < 20) |>
    post_request()
  feats <- items$features
  info  <- data.frame(i = seq_along(feats),
                      tile = sapply(feats, function(f) substr(f$id, 5, 10)),
                      cc = sapply(feats, function(f) f$properties$`eo:cloud_cover`))
  keep <- info |> group_by(tile) |> arrange(cc, .by_group = TRUE) |> slice_head(n = 2) |> pull(i)
  feats <- feats[keep]
  cat(length(feats), "Sentinel-2 scenes used:\n"); print(sapply(feats, function(f) f$id))

  bands <- c(blue = "blue", green = "green", red = "red", nir = "nir", swir = "swir16")
  scenes <- lapply(feats, function(f) tryCatch({
    cat("  reading", f$id, format(Sys.time(), "%H:%M"), "\n")
    b <- rast(lapply(bands, function(k)
      warp_to_grid(f$assets[[k]]$href, sprintf("data/raw/s2_30m/%s_%s.tif", f$id, k), "average", 0)))
    scl <- warp_to_grid(f$assets$scl$href, sprintf("data/raw/s2_30m/%s_scl.tif", f$id), "mode", 0)
    if (global(b[[1]], "max", na.rm = TRUE)[1, 1] > 2) b <- b * 0.0001 - 0.1   # scale only if GDAL has not already
    mask(b, scl %in% c(0, 1, 2, 3, 8, 9, 10, 11), maskvalues = TRUE)   # clouds, shadow, no data
  }, error = function(e) { message("  skipped ", f$id, ": ", conditionMessage(e)); NULL }))
  scenes <- Filter(function(x) !is.null(x) && global(!is.na(x[[1]]), "sum")[1, 1] > 0, scenes)
  stopifnot(length(scenes) >= 1)
  cat(length(scenes), "scenes with data composited\n")
  # per-band median over scenes
  s2 <- rast(lapply(1:5, function(i) app(rast(lapply(scenes, `[[`, i)), median, na.rm = TRUE)))
  names(s2) <- names(bands)
  writeRaster(s2, "data/s2_composite.tif", overwrite = TRUE)
}
s2 <- align(rast("data/s2_composite.tif"))
ndvi  <- (s2$nir - s2$red) / (s2$nir + s2$red);                        names(ndvi)  <- "ndvi"
bsi   <- ((s2$swir + s2$red) - (s2$nir + s2$blue)) /
         ((s2$swir + s2$red) + (s2$nir + s2$blue));                     names(bsi)   <- "bsi"
mndwi <- (s2$green - s2$swir) / (s2$green + s2$swir);                  names(mndwi) <- "mndwi"

# 2c. Land cover (ESA WorldCover 2021, 10 m -> 30 m majority) ----------------
cat("Land cover", format(Sys.time(), "%H:%M"), "\n")
lc_parts <- lapply(c("N24E090", "N24E093"), function(t) tryCatch(
  warp_to_grid(sprintf("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_%s_Map.tif", t),
               sprintf("data/raw/worldcover_%s.tif", t), "mode", 0),
  error = function(e) NULL))
lc_parts <- Filter(Negate(is.null), lc_parts)
lc_parts <- lapply(lc_parts, align, method = "near")
lc <- if (length(lc_parts) > 1) cover(lc_parts[[1]], lc_parts[[2]]) else lc_parts[[1]]
names(lc) <- "lc"

# 2d. Flood frequency (JRC Global Surface Water occurrence) -----------------
cat("Water occurrence", format(Sys.time(), "%H:%M"), "\n")
USE_JRC <- TRUE      # set to FALSE if this download is too slow on your connection
occ <- NULL
if (USE_JRC) {
  occ <- tryCatch(warp_to_grid("https://storage.googleapis.com/global-surface-water/downloads2021/occurrence/occurrence_90E_30Nv1_4_2021.tif",
                               "data/water_occurrence.tif", "average"),
                  error = function(e) { message("JRC occurrence not available: ", conditionMessage(e)); NULL })
  if (!is.null(occ)) { occ <- classify(align(occ), cbind(NA, 0)); names(occ) <- "water_occ" }
}

# 2e. Access: distance to road (OpenStreetMap) -----------------------------
cat("Roads", format(Sys.time(), "%H:%M"), "\n")
if (!file.exists("data/roads.gpkg")) {
  library(osmdata)
  rd <- opq(st_bbox(aoiG), timeout = 300) |>
    add_osm_feature("highway", c("motorway", "trunk", "primary", "secondary", "tertiary",
                                 "unclassified", "residential", "track", "service")) |>
    osmdata_sf()
  st_write(st_transform(rd$osm_lines, 32646)[, "highway"], "data/roads.gpkg", delete_dsn = TRUE, quiet = TRUE)
}
roads <- st_read("data/roads.gpkg", quiet = TRUE)
droad <- distance(rasterize(vect(roads), tmpl)); names(droad) <- "dist_road"

# 2f. Stack, mask non-soil areas, save ------------------------------------
cov <- c(dem, slope, tpi, twi, hand, dstream, ndvi, bsi, mndwi, lc, droad)
if (!is.null(occ)) cov <- c(cov, occ)
cov <- mask(cov, vect(ws))

soil <- !(cov$lc %in% c(50, 80))                     # drop built-up and permanent water
if (!is.null(occ)) soil <- soil & (cov$water_occ < 50)
soil <- soil & (cov$mndwi < 0.2)                     # drop open water seen in the composite
cov_soil <- mask(cov, soil, maskvalues = FALSE)

writeRaster(cov, "data/covariates_all.tif", overwrite = TRUE)
writeRaster(cov_soil, "data/covariates_soil.tif", overwrite = TRUE)
writeRaster(streams, "data/streams_30m.tif", overwrite = TRUE)

# 2g. Summary table and figure --------------------------------------------
sm <- global(cov_soil, c("min", "mean", "max"), na.rm = TRUE)
sm$covariate <- rownames(sm)
write.csv(sm, "outputs/covariate_summary.csv", row.names = FALSE)
cat(sprintf("Soil area kept: %.0f ha of %.0f ha\n",
            global(!is.na(cov_soil$elev), "sum")[1, 1] * RES^2 / 1e4,
            as.numeric(st_area(ws)) / 1e4))

png("outputs/figures/02_covariates.png", width = 2400, height = 1800, res = 200)
plot(cov_soil[[setdiff(names(cov_soil), "lc")]], nc = 4, mar = c(1, 1, 2, 4), axes = FALSE,
     col = hcl.colors(50, "viridis"))
dev.off()
