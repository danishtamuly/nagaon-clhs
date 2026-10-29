# =============================================================================
# Step 4: Figures for the report and presentation
# =============================================================================
library(sf); library(terra); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
sf_use_s2(FALSE)
load("outputs/clhs_results.RData")

ws    <- st_read("data/watershed.gpkg", quiet = TRUE)
roads <- st_read("data/roads.gpkg", quiet = TRUE) |> st_intersection(st_buffer(ws, 200))
roads <- roads |> filter(highway %in% c("trunk", "primary", "secondary", "tertiary", "unclassified", "residential"))
dem   <- rast("data/dem.tif") |> mask(vect(ws))
cov   <- rast("data/covariates_all.tif")
nagaon <- st_read("data/nagaon_district.gpkg", quiet = TRUE) |> st_transform(32646)

th <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_line(colour = "grey92"), plot.title = element_text(face = "bold"))
col_site <- "#C8384E"; col_back <- "#B7791F"

# 4a. Location map ------------------------------------------------------------
p_loc <- ggplot() +
  geom_sf(data = nagaon, fill = "grey95", colour = "grey40") +
  geom_sf(data = ws, fill = col_site, colour = col_site, alpha = 0.6) +
  labs(title = "Study watershed in Nagaon district") + th
ggsave("outputs/figures/04a_location.png", p_loc, width = 6, height = 5, dpi = 220)

# 4b. Site map on HAND with hillshade ------------------------------------------
hs <- shade(terrain(dem, "slope", unit = "radians"), terrain(dem, "aspect", unit = "radians"), 40, 315)
hand_df <- as.data.frame(mask(cov$hand, vect(ws)), xy = TRUE, na.rm = TRUE)
hs_df   <- as.data.frame(hs, xy = TRUE, na.rm = TRUE); names(hs_df)[3] <- "hs"
prof <- out |> filter(role == "profile"); back <- out |> filter(role == "replacement")
p_map <- ggplot() +
  geom_raster(data = hand_df, aes(x, y, fill = pmin(hand, 15))) +
  scale_fill_gradientn(colours = c("#0B2A5B", "#2E7D8C", "#8FC1A9", "#F1E8C8"), name = "HAND (m)") +
  geom_raster(data = hs_df, aes(x, y, alpha = hs), fill = "black", show.legend = FALSE) +
  scale_alpha(range = c(0.25, 0)) +
  geom_sf(data = roads, colour = "grey25", linewidth = 0.2) +
  geom_sf(data = ws, fill = NA, colour = "black", linewidth = 0.6) +
  geom_point(data = back, aes(x, y), shape = 21, fill = "white", colour = col_back, size = 1.6, stroke = 0.8) +
  geom_point(data = prof, aes(x, y), shape = 21, fill = col_site, colour = "white", size = 3, stroke = 0.8) +
  geom_text(data = prof, aes(x, y, label = sub(".*-P", "", site_id)), size = 2.6, nudge_y = 450, fontface = "bold") +
  coord_sf(crs = 32646) +
  labs(title = "30 cLHS soil profile sites", subtitle = "Red: profile sites   Hollow: replacement sites   Lines: roads", x = NULL, y = NULL) + th
ggsave("outputs/figures/04b_site_map.png", p_map, width = 8, height = 8, dpi = 220)

# 4c. Coverage of each covariate: watershed vs sites ----------------------------
labs_v <- c(elev = "Elevation (m)", slope = "Slope (deg)", tpi = "TPI (m)", twi = "Wetness index",
            hand = "HAND (m)", dist_stream = "Distance to stream (m)", ndvi = "NDVI (rabi)",
            bsi = "Bare soil index", mndwi = "MNDWI", water_occ = "Water occurrence (%)")
long_pop  <- pool[, cont]  |> pivot_longer(everything(), names_to = "v")
long_site <- sites[, cont] |> pivot_longer(everything(), names_to = "v")
p_cov <- ggplot(long_pop, aes(value)) +
  geom_histogram(bins = 30, fill = "grey80", colour = NA) +
  geom_rug(data = long_site, colour = col_site, linewidth = 0.7, length = unit(0.12, "npc")) +
  facet_wrap(~ v, scales = "free", ncol = 5, labeller = as_labeller(labs_v)) +
  labs(title = "Do the 30 sites span every covariate?",
       subtitle = "Grey: whole watershed (candidate points). Red ticks: cLHS sites.", x = NULL, y = NULL) + th
ggsave("outputs/figures/04c_coverage.png", p_cov, width = 12, height = 5, dpi = 220)

# 4d. Land cover share ---------------------------------------------------------
p_lc <- ggplot(lc_tab, aes(reorder(lc_name, pct), pct, fill = set)) +
  geom_col(position = "dodge", width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("Watershed" = "grey70", "cLHS sites" = col_site), name = NULL) +
  labs(title = "Land cover: watershed vs sites", x = NULL, y = "% of points") + th
ggsave("outputs/figures/04d_landcover.png", p_lc, width = 6, height = 4, dpi = 220)

# 4e. cLHS vs random designs (Monte Carlo) -------------------------------------
p_mc <- ggplot(data.frame(cov = 100 * mc), aes(cov)) +
  geom_histogram(bins = 25, fill = "grey75") +
  geom_vline(xintercept = 100 * mean(diag$strata_filled), colour = col_site, linewidth = 1.2) +
  annotate("text", x = 100 * mean(diag$strata_filled), y = Inf, vjust = 1.5, hjust = 1.05,
           label = sprintf("cLHS: %.0f%%", 100 * mean(diag$strata_filled)), colour = col_site, fontface = "bold") +
  labs(title = "cLHS vs 500 random designs of 30 sites",
       x = "Average % of covariate strata filled", y = "Number of random designs") + th
ggsave("outputs/figures/04e_clhs_vs_random.png", p_mc, width = 7, height = 4, dpi = 220)

cat("Figures written to outputs/figures\n")
