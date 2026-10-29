# =============================================================================
# Step 3: Select 30 soil profile sites by conditioned Latin hypercube sampling
# Method: Minasny & McBratney (2006), clhs package (Roudier et al.),
#         candidates limited to land within 1 km of a road; 3 replacement sites per profile.
# =============================================================================
library(sf); library(terra); library(dplyr); library(tidyr); library(ggplot2); library(clhs)
sf_use_s2(FALSE)

N_PROFILES <- 30
ITER       <- 50000
ACCESS_M   <- 1000     # only consider land within 1 km of a road (field access)
MIN_KM     <- 1.5        # minimum spacing between profiles
N_BACKUP   <- 3
BACKUP_KM  <- 1.0
PREFIX     <- "NGW"      # Nagaon Watershed
set.seed(2026)

cov <- rast("data/covariates_soil.tif")
ws  <- st_read("data/watershed.gpkg", quiet = TRUE)
cont <- intersect(c("elev", "slope", "tpi", "twi", "hand", "dist_stream",
                    "ndvi", "bsi", "mndwi"), names(cov))
# water_occ (flood frequency) is zero on ~90% of the land, so it cannot be cut into
# 30 equal strata; it is used to mask permanent water in step 2 instead.

# 3a. Candidate points: regular grid sample of the soil area -----------------
pool <- spatSample(cov, size = 12000, method = "regular", xy = TRUE, na.rm = TRUE) |>
  as.data.frame() |> tidyr::drop_na()
pool$lc <- factor(pool$lc)
n_all <- nrow(pool)
pool <- pool[pool$dist_road <= ACCESS_M, ]
cat(nrow(pool), "candidate points within", ACCESS_M, "m of a road (of", n_all, "in the soil area)\n")
x <- pool[, c(cont, "lc")]

# 3b. Run cLHS with several seeds; keep the best design that respects spacing --
spacing_ok <- function(idx) {
  d <- as.matrix(dist(pool[idx, c("x", "y")])) / 1000
  sum(d[upper.tri(d)] < MIN_KM)
}
best <- NULL
for (s in 1:5) {
  set.seed(2026 + s)
  r <- clhs(x, size = N_PROFILES, iter = ITER, simple = FALSE, progress = FALSE)
  v <- spacing_ok(r$index_samples)
  o <- tail(r$obj, 1)
  cat(sprintf("run %d: objective %.2f, pairs closer than %.1f km: %d\n", s, o, MIN_KM, v))
  if (is.null(best) || (v < best$v) || (v == best$v && o < best$o)) best <- list(res = r, v = v, o = o, seed = 2026 + s)
}
res <- best$res
idx <- res$index_samples
cat("Chosen run seed", best$seed, "\n")

# Spacing fix: if two sites are closer than MIN_KM, move the second one to the
# most similar candidate (same land cover, nearest in covariate space) that is
# at least MIN_KM from every other site.
z <- scale(pool[, cont])
moved <- 0
repeat {
  d <- as.matrix(dist(pool[idx, c("x", "y")])) / 1000; diag(d) <- Inf
  bad <- which(d < MIN_KM, arr.ind = TRUE)
  if (!nrow(bad)) break
  j <- max(bad[1, ]); s0 <- idx[j]
  others <- pool[idx[-j], c("x", "y")]
  far <- sapply(seq_len(nrow(pool)), function(k) min(sqrt((others$x - pool$x[k])^2 + (others$y - pool$y[k])^2))) / 1000
  ok <- which(far >= MIN_KM & pool$lc == pool$lc[s0] & !(seq_len(nrow(pool)) %in% idx))
  if (!length(ok)) ok <- which(far >= MIN_KM & !(seq_len(nrow(pool)) %in% idx))
  dz <- rowSums((z[ok, , drop = FALSE] - matrix(z[s0, ], length(ok), ncol(z), byrow = TRUE))^2)
  idx[j] <- ok[which.min(dz)]; moved <- moved + 1
}
cat("Sites moved to keep", MIN_KM, "km spacing:", moved, "\n")

png("outputs/figures/03a_clhs_objective.png", width = 1600, height = 900, res = 200)
plot(res, modes = "obj")
dev.off()

sites <- pool[idx, ] |>
  mutate(site_id = sprintf("%s-P%02d", PREFIX, row_number()), role = "profile")

# 3c. Replacement sites: same land cover, within 1 km, nearest in covariate space
used <- idx
backs <- list()
for (i in seq_along(idx)) {
  s  <- idx[i]
  dk <- sqrt((pool$x - pool$x[s])^2 + (pool$y - pool$y[s])^2) / 1000
  ok <- which(dk <= BACKUP_KM & dk > 0.1 & pool$lc == pool$lc[s] & !(seq_len(nrow(pool)) %in% used))
  if (!length(ok)) next
  dz <- sqrt(rowMeans((z[ok, , drop = FALSE] - matrix(z[s, ], length(ok), length(cont), byrow = TRUE))^2))
  pick <- ok[order(dz)][seq_len(min(N_BACKUP, length(ok)))]
  used <- c(used, pick)
  backs[[i]] <- pool[pick, ] |>
    mutate(site_id = paste0(sites$site_id[i], letters[seq_along(pick)]),
           role = "replacement", parent_id = sites$site_id[i],
           dist_km = round(dk[pick], 2), cov_dist_sd = round(sort(dz)[seq_along(pick)], 2))
}
out <- bind_rows(sites, bind_rows(backs))

# 3d. Coordinates and exports ---------------------------------------------
pts <- st_as_sf(out, coords = c("x", "y"), crs = 32646, remove = FALSE)
ll  <- st_coordinates(st_transform(pts, 4326))
out$lon <- round(ll[, 1], 6); out$lat <- round(ll[, 2], 6)
lc_names <- c(`10` = "Tree cover", `20` = "Shrubland", `30` = "Grassland", `40` = "Cropland",
              `50` = "Built-up", `60` = "Bare", `80` = "Water", `90` = "Herbaceous wetland")
out$lc_name <- unname(lc_names[as.character(out$lc)])
out <- out |> relocate(site_id, role, any_of("parent_id"), lat, lon, lc_name)
write.csv(out, "outputs/field/nagaon_ws_clhs_sites.csv", row.names = FALSE)

pts <- st_as_sf(out, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
st_write(pts, "outputs/field/nagaon_ws_clhs_sites.gpkg", delete_dsn = TRUE, quiet = TRUE)
kml_pm <- sprintf(
  "<Placemark><name>%s</name><styleUrl>#%s</styleUrl><description>%s; %s; elevation %.0f m; HAND %.1f m</description><Point><coordinates>%.6f,%.6f,0</coordinates></Point></Placemark>",
  out$site_id, out$role, out$role, out$lc_name, out$elev, out$hand, out$lon, out$lat)
writeLines(c('<?xml version="1.0" encoding="UTF-8"?>',
  '<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>Nagaon watershed cLHS sites</name>',
  '<Style id="profile"><IconStyle><color>ff4e38c8</color><scale>1.1</scale></IconStyle></Style>',
  '<Style id="replacement"><IconStyle><color>ff1f79b7</color><scale>0.7</scale></IconStyle></Style>',
  kml_pm, '</Document></kml>'), "outputs/field/nagaon_ws_clhs_sites.kml")
st_write(transmute(pts, name = site_id, desc = paste(role, lc_name)),
         "outputs/field/nagaon_ws_clhs_sites.gpx", driver = "GPX", layer = "waypoints",
         delete_dsn = TRUE, quiet = TRUE)

# 3e. How well do the 30 sites cover the watershed? --------------------------
strata_cover <- function(sx, px, n) {
  br <- quantile(px, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)
  length(unique(findInterval(sx, br[-c(1, n + 1)]))) / n
}
ks_d <- function(sx, px) suppressWarnings(ks.test(sx, px)$statistic)
diag <- data.frame(
  covariate = cont,
  strata_filled = sapply(cont, function(v) strata_cover(sites[[v]], pool[[v]], N_PROFILES)),
  ks_D = sapply(cont, function(v) ks_d(sites[[v]], pool[[v]])),
  pop_median = sapply(cont, function(v) median(pool[[v]])),
  sample_median = sapply(cont, function(v) median(sites[[v]])))
write.csv(diag, "outputs/clhs_diagnostics.csv", row.names = FALSE)
print(diag)

# Monte Carlo comparison: 500 simple random designs of the same size
mc <- replicate(500, {
  r <- sample(nrow(pool), N_PROFILES)
  mean(sapply(cont, function(v) strata_cover(pool[[v]][r], pool[[v]], N_PROFILES)))
})
mc_summary <- data.frame(design = c("cLHS (this study)", "Random: mean of 500", "Random: best of 500"),
                         mean_strata_filled = c(mean(diag$strata_filled), mean(mc), max(mc)))
write.csv(mc_summary, "outputs/clhs_vs_random.csv", row.names = FALSE)
print(mc_summary)

lc_tab <- bind_rows(
  data.frame(set = "Watershed", lc = pool$lc),
  data.frame(set = "cLHS sites", lc = sites$lc)) |>
  count(set, lc) |> group_by(set) |> mutate(pct = 100 * n / sum(n)) |> ungroup() |>
  mutate(lc_name = unname(lc_names[as.character(lc)]))
write.csv(lc_tab, "outputs/landcover_share.csv", row.names = FALSE)

save(pool, sites, out, diag, mc, mc_summary, lc_tab, cont, file = "outputs/clhs_results.RData")
cat("Done:", nrow(sites), "profiles,", sum(out$role == "replacement"), "replacement sites\n")
