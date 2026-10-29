# =============================================================================
# Step 1: Choose the study watershed in Nagaon
# Data: HydroBASINS levels 9-10 (Lehner & Grill 2013), GADM 4.1 district
#       boundary, Copernicus GLO-30 DEM.
# Output: outputs/watershed_candidates.csv, outputs/figures/01_watershed_candidates.png
# =============================================================================
library(sf); library(terra); library(dplyr); library(ggplot2)
sf_use_s2(FALSE)
options(timeout = 3600)

# 1a. Nagaon district boundary (GADM 4.1) -------------------------------------
ind <- geodata::gadm("IND", level = 2, path = "data/raw")
nagaon <- st_as_sf(ind) |> filter(NAME_1 == "Assam", grepl("Nagaon|Nowgong", NAME_2))
stopifnot(nrow(nagaon) >= 1)
nagaon <- nagaon |> summarise()
st_write(nagaon, "data/nagaon_district.gpkg", delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("Nagaon district: %.0f km2\n", as.numeric(st_area(st_transform(nagaon, 32646))) / 1e6))

# 1b. HydroBASINS sub-watersheds (Asia, levels 9 and 10) -----------------------
zip <- "data/raw/hybas_as_lev01-12_v1c.zip"
if (!file.exists(zip))
  download.file("https://data.hydrosheds.org/file/hydrobasins/standard/hybas_as_lev01-12_v1c.zip",
                zip, mode = "wb")
shp <- grep("lev(09|10)_v1c\\.(shp|shx|dbf|prj)$", unzip(zip, list = TRUE)$Name, value = TRUE)
unzip(zip, files = shp, exdir = "data/raw/hybas")
read_lev <- function(l) {
  f <- list.files("data/raw/hybas", sprintf("lev%s_v1c\\.shp$", l), full.names = TRUE, recursive = TRUE)
  st_read(f, quiet = TRUE) |> st_crop(st_bbox(st_buffer(nagaon, 0.05))) |> mutate(level = as.integer(l))
}
hb <- bind_rows(read_lev("09"), read_lev("10"))

# keep sub-watersheds 100-500 km2 lying >= 90% inside Nagaon
inter <- st_intersection(hb, nagaon)
inter$in_km2 <- as.numeric(st_area(st_transform(inter, 32646))) / 1e6
cand <- hb |>
  left_join(st_drop_geometry(inter)[, c("HYBAS_ID", "level", "in_km2")], by = c("HYBAS_ID", "level")) |>
  mutate(share_in = in_km2 / SUB_AREA) |>
  filter(share_in >= 0.9, SUB_AREA >= 100, SUB_AREA <= 500)
cat(nrow(cand), "candidate sub-watersheds\n")

# 1c. Terrain variety in each candidate (Copernicus GLO-30 DEM) ------------------
tiles <- c("N25_00_E092_00", "N25_00_E093_00", "N26_00_E092_00", "N26_00_E093_00")
urls  <- sprintf("/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_%s_DEM/Copernicus_DSM_COG_10_%s_DEM.tif", tiles, tiles)
box   <- vect(st_buffer(nagaon, 0.03))
dem   <- merge(sprc(lapply(urls, function(u) crop(rast(u), box))))
writeRaster(dem, "data/raw/dem_nagaon_cop30.tif", overwrite = TRUE)

ex <- terra::extract(dem, vect(cand))
names(ex) <- c("ID", "z")
st <- ex |> group_by(ID) |>
  summarise(elev_min = min(z, na.rm = TRUE), elev_p95 = quantile(z, 0.95, na.rm = TRUE),
            elev_sd = sd(z, na.rm = TRUE), low_share = mean(z < min(z, na.rm = TRUE) + 5, na.rm = TRUE))
cand <- cand |> mutate(ID = row_number()) |> left_join(st, by = "ID") |>
  mutate(relief = elev_p95 - elev_min) |>
  arrange(desc(relief)) |> mutate(rank = row_number())

tab <- st_drop_geometry(cand) |>
  transmute(rank, HYBAS_ID, level, area_km2 = round(SUB_AREA), elev_min = round(elev_min),
            elev_p95 = round(elev_p95), relief = round(relief), elev_sd = round(elev_sd, 1),
            lowland_pct = round(100 * low_share))
print(tab)
write.csv(tab, "outputs/watershed_candidates.csv", row.names = FALSE)
st_write(cand, "data/watershed_candidates.gpkg", delete_dsn = TRUE, quiet = TRUE)

p <- ggplot() +
  geom_sf(data = nagaon, fill = "grey96", colour = "grey30") +
  geom_sf(data = cand, aes(fill = relief), colour = "white", linewidth = 0.3, alpha = 0.9) +
  geom_sf_text(data = cand, aes(label = rank), size = 3, fontface = "bold") +
  scale_fill_viridis_c(name = "Relief (m)") +
  labs(title = "Candidate sub-watersheds in Nagaon",
       subtitle = "HydroBASINS levels 9-10, 100-500 km2, >=90% inside the district; numbered by relief",
       x = NULL, y = NULL) +
  theme_minimal()
ggsave("outputs/figures/01_watershed_candidates.png", p, width = 8, height = 6.5, dpi = 200)
cat("Step 1 done. Open outputs/figures/01_watershed_candidates.png\n")
