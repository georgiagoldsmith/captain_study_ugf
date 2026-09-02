# Preparing the input layers for CAPTAIN
# G Goldsmtih 
# August 2026
###############################################################################

# Cost, disturbance and protected areas are taken from the layers prioritizr
# wrote to outputs/layers_3km, so both models are fed the same data.
# RUN scripts/prioritizr_ugf_prioritization.R FIRST -- this script stops if
# those layers are missing.

library(here)
library(terra)
library(sf)

OUT <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected"
dir.create(file.path(OUT, "present_habitat_suitability"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "environmental_layers"),        recursive = TRUE, showWarnings = FALSE)
ENV <- file.path(OUT, "environmental_layers")
SDM <- file.path(OUT, "present_habitat_suitability")

ugf   <- vect(st_transform(st_read(here("data/UGF_gp.shp","UGF_gp.shp"), quiet = TRUE), 3857))
# The grid is taken from an existing CAPTAIN export rather than rebuilt, so
# results stay comparable with every run done before this one.
grid  <- rast("/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data/environmental_layers/costs.tif")
grid  <- ifel(!is.na(grid), 0, grid)
cat(sprintf("3km grid: %d x %d, res %.0f m\n", nrow(grid), ncol(grid), res(grid)[1]))

# Layers built by the prioritizr script. They are on the same 3km grid, so they
# are read as-is; compareGeom() below stops the run if that ever stops being true.
LAYERS <- here("outputs/layers_3km")
if (!dir.exists(LAYERS))
  stop("outputs/layers_3km not found -- run scripts/prioritizr_ugf_prioritization.R first")

# --- cost ---
# Cost is cocoa suitability, already 0-1, so this only aligns it to the grid.
# Cells with no value are places the cocoa model made no prediction, mostly
# urban, water and cropland. They are filled with 0 rather than dropped.
cost <- rast(file.path(LAYERS, "cocoa.tif"))
stopifnot(compareGeom(cost, grid, stopOnError = FALSE))
inside <- !is.na(rasterize(ugf, grid, field = 1, touches = TRUE))
n_absent <- global(is.na(cost) & inside, "sum", na.rm = TRUE)[[1]]
cost <- ifel(is.na(cost) & inside, 0, cost)      # zero suitability, zero cost
cost <- mask(cost, ugf, touches = TRUE)
cat(sprintf("cost: filled %d absent cells with 0 (zero cocoa suitability)\n", n_absent))
writeRaster(cost, file.path(ENV, "costs.tif"), overwrite = TRUE)

# --- disturbance ---
# gHM is a measured, continuous score of human pressure, so it is used directly
# as disturbance. Urban and water are set to exactly 1.0 -- "not selectable".
lccs <- mask(project(rast(here("data/LC/lccs_ugf.tif")), grid, method = "near"), ugf, touches = TRUE)
ghm  <- rast(file.path(LAYERS, "ghm.tif"))
stopifnot(compareGeom(ghm, grid, stopOnError = FALSE))
ghm  <- mask(ghm, ugf, touches = TRUE)
exclude_classes <- c(190, 210)                       # urban, water
disturbance <- ifel(lccs %in% exclude_classes, 1, ghm)
# adjust disturbance levels of 1 to 0.999 so they are still selectable
disturbance <- ifel(!(lccs %in% exclude_classes) & disturbance >= 1, 0.999, disturbance)
writeRaster(disturbance, file.path(ENV, "disturbance.tif"), overwrite = TRUE)
# gHM on its own, without the urban/water override, kept for reference
writeRaster(ghm, file.path(ENV, "disturbance_ghm_raw.tif"), overwrite = TRUE)

# --- protected areas ---
# Rasterized protected areas, where cells count as protected only if its center
# falls inside a protected area boundary.
pa <- rast(file.path(LAYERS, "protected_areas.tif"))
stopifnot(compareGeom(pa, grid, stopOnError = FALSE))
pa <- mask(pa, ugf, touches = TRUE)
writeRaster(pa, file.path(ENV, "protected_areas.tif"), overwrite = TRUE)

# --- species habitat ---
# One 0/1 raster per species, read from the same folder prioritizr uses so both
# models get identical habitat. 
fs <- list.files(here("data/birds_aoh/clipped"), pattern = "tif$", full.names = TRUE)
for (f in fs) {
  r <- project(rast(f), grid, method = "near")
  # These files store no-data in an integer type that cannot hold it, so
  # reprojecting turns empty cells into a huge number. Real habitat is always
  # exactly 1, so treat anything else as no-data.
  r <- ifel(r == 1, 1, NA)
  r <- mask(r, ugf, touches = TRUE)
  # Two of the habitat files have an "_R" suffix that species_traits.csv does
  # not. CAPTAIN pairs them up by filename, so drop the suffix here.
  writeRaster(r, file.path(SDM, sub("_R\\.tif$", ".tif", basename(f))), overwrite = TRUE)
}

cat(sprintf("\nwrote %d species + 4 env layers to\n  %s\n", length(fs), OUT))
for (f in c("costs.tif","disturbance.tif","protected_areas.tif")) {
  r <- rast(file.path(ENV,f)); v <- values(r); v <- v[!is.na(v) & !is.nan(v)]
  cat(sprintf("  %-22s %dx%d  valid %7d  min %.3f  max %.3f  mean %.3f\n",
              f, nrow(r), ncol(r), length(v), min(v), max(v), mean(v)))
}
