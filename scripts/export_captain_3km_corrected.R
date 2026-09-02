# Export CAPTAIN inputs on the 3km grid, built from the CORRECT sources.
#
# The original 3km export targeted `cocoa_ugf` -- a coarse 9.6km raster that is
# NOT cocoa suitability, and the source of the costs.tif bug (prioritzr.R line
# 208 writes cocoa_ugf_norm, never assigned in that script, instead of line
# 161's cocoa_suitability_ugf_norm). This rebuilds the same 3km grid from the
# right sources: real cocoa suitability for cost, gHM for disturbance, urban and
# water only for hard exclusion.
#
# 3km keeps the dispersal cutoff of 3 cells = 9km, so all 18 species retain
# their AVONET Hand-Wing dispersal rates. At 1km that combination needs 22 GB
# per runner and does not fit in 24 GiB.
#
#   Rscript scripts/export_captain_1km.R

library(here); library(terra); library(sf)

OUT <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected"
dir.create(file.path(OUT, "present_habitat_suitability"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "environmental_layers"),        recursive = TRUE, showWarnings = FALSE)
ENV <- file.path(OUT, "environmental_layers")
SDM <- file.path(OUT, "present_habitat_suitability")

ugf   <- vect(st_transform(st_read(here("data/UGF_gp.shp","UGF_gp.shp"), quiet = TRUE), 3857))
# the 3km grid CAPTAIN uses, taken from the existing export so every previous
# run stays directly comparable
grid  <- rast("/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data/environmental_layers/costs.tif")
grid  <- ifel(!is.na(grid), 0, grid)
cat(sprintf("3km grid: %d x %d, res %.0f m\n", nrow(grid), ncol(grid), res(grid)[1]))

# --- cost: already normalised 0-1 at 963m, just align and mask -------------
# Cells with NO cost value are where the cocoa model made no prediction -- mostly
# urban / water / cropland. Cost here IS cocoa suitability, and you cannot grow
# cocoa on water or tarmac, so their suitability is genuinely 0. They are kept
# out of selection by the disturbance layer, not by an inflated cost (the same
# separation prioritzr.R line 206 describes).
cost <- resample(rast(here("data/cocoa suitability","cost.tif")), grid, method = "bilinear")
inside <- !is.na(rasterize(ugf, grid, field = 1, touches = TRUE))
n_absent <- global(is.na(cost) & inside, "sum", na.rm = TRUE)[[1]]
cost <- ifel(is.na(cost) & inside, 0, cost)      # zero suitability, zero cost
cost <- mask(cost, ugf, touches = TRUE)
cat(sprintf("cost: filled %d absent cells with 0 (zero cocoa suitability)\n", n_absent))
writeRaster(cost, file.path(ENV, "costs.tif"), overwrite = TRUE)

# --- disturbance: gHM everywhere, urban/water as the exclusion sentinel ----
# gHM is a continuous, measured human-modification index, so it replaces the
# LCCS tiering entirely rather than only filling the cropland/bare classes.
# LCCS is still consulted, but only to mark the two classes that cannot be
# protected at all (urban, water) with the 1.0 sentinel CAPTAIN's action mask
# keys off. Cropland and bare are no longer excluded -- they carry their gHM
# value and stay selectable.
lccs <- mask(project(rast(here("data/LC/lccs_ugf.tif")), grid, method = "near"), ugf, touches = TRUE)
ghm  <- mask(project(rast(here("data/gHM/gHM.tif")), grid, method = "average"), ugf, touches = TRUE)
exclude_classes <- c(190, 210)                       # urban, water
disturbance <- ifel(lccs %in% exclude_classes, 1, ghm)
# guard: a gHM reading of exactly 1.0 outside those classes would be silently
# excluded by the >= 1.0 action mask
disturbance <- ifel(!(lccs %in% exclude_classes) & disturbance >= 1, 0.999, disturbance)
writeRaster(disturbance, file.path(ENV, "disturbance.tif"), overwrite = TRUE)
# raw gHM without the sentinel, kept for reference
writeRaster(ghm, file.path(ENV, "disturbance_ghm_raw.tif"), overwrite = TRUE)

# --- protected areas: rasterized from the shapefile ------------------------
# centre-based (no touches=TRUE), per prioritzr.R lines 38-42
pa <- rasterize(vect(st_read(here("data/protected areas","protected_areas_ugf.shp"), quiet = TRUE)),
                grid, field = 1, background = 0)
pa <- mask(pa, ugf, touches = TRUE)
writeRaster(pa, file.path(ENV, "protected_areas.tif"), overwrite = TRUE)

# --- species AOH, normalised 0-1 ------------------------------------------
fs <- list.files(here("data/birds_aoh/birds_aoh"), pattern = "tif$", full.names = TRUE)
for (f in fs) {
  r <- project(rast(f), grid, method = "near")
  r <- ifel(r == 1, 1, NA)                 # guard the INT4U overflow, prioritzr.R lines 145-149
  r <- mask(r, ugf, touches = TRUE)
  # AOH filenames carry an "_R" (resident) suffix for two species; species_traits.csv
  # does not, and CAPTAIN matches the two by name. Strip it here.
  writeRaster(r, file.path(SDM, sub("_R\\.tif$", ".tif", basename(f))), overwrite = TRUE)
}

cat(sprintf("\nwrote %d species + 4 env layers to\n  %s\n", length(fs), OUT))
for (f in c("costs.tif","disturbance.tif","protected_areas.tif")) {
  r <- rast(file.path(ENV,f)); v <- values(r); v <- v[!is.na(v) & !is.nan(v)]
  cat(sprintf("  %-22s %dx%d  valid %7d  min %.3f  max %.3f  mean %.3f\n",
              f, nrow(r), ncol(r), length(v), min(v), max(v), mean(v)))
}
