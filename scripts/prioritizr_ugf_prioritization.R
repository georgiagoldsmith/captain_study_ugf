# Conservation prioritisation for the Upper Guinean Forest, on a 3km grid.
#
#   Rscript scripts/prioritizr_ugf_prioritization.R
#
# Chooses roughly 30% of the study area to protect: the cells that best cover
# 18 bird species' habitat while steering away from land that is already
# degraded or valuable for growing cocoa. Existing protected areas are locked
# in, urban land is locked out. Solved as an integer program (prioritizr +
# Gurobi).
#
# Two ideas do the work, and they are deliberately different mechanisms:
#
#   HABITAT QUALITY, applied as a discount. Human pressure (gHM, the Global
#   Human Modification index: 0 = untouched, 1 = fully modified) multiplies
#   each species' habitat layer. A degraded cell therefore still costs a full
#   unit of budget but contributes less toward that species' target. This can
#   never price a cell out of the solution -- a disturbed cell that happens to
#   be the only habitat for a critically endangered bird still gets picked.
#
#   OPPORTUNITY COST, applied as a penalty. Cocoa suitability is charged
#   against the objective, so choosing land that is good for cocoa actively
#   costs something. Unlike the discount, this CAN price a cell out.
#
# One formulation, solved twice: once with cocoa free, to measure the baseline
# the penalty is scaled against, then once for real. COCOA_SHARE sets how hard
# cocoa is charged for.
#
# INPUTS   data/                    UGF boundary, cocoa suitability, gHM,
#                                   protected areas, urban land cover, 18
#                                   species habitat rasters, species traits
#          captain/disturbance.tif  read for its grid geometry only
#
# OUTPUTS  outputs/prioritizr_p3a_ghm_discount_3km.tif  the solution
#          outputs/p3a_ghm_discount_3km_coverage.csv    coverage by species
#          outputs/solution_3a_ghm_discount_3km.png     map of the solution
library(here)
library(terra)
library(sf)
library(prioritizr)

t_start <- Sys.time()
stamp <- function(msg) cat(sprintf("[%5.1fs] %s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs")), msg))

# ---- grid ----
# Everything is built on the same 3km grid the CAPTAIN model outputs use, so
# the two can be compared cell for cell. Only the geometry of disturbance.tif
# is used, never its values. Note the cells are not square: 3001.4 x 2993.9 m.
template <- rast(here("captain/disturbance.tif"))
ugf_boundary <- st_read(here("data/UGF_gp.shp", "UGF_gp.shp"), quiet = TRUE)
ugf_vect_3857 <- vect(st_transform(ugf_boundary, 3857))

# ---- cocoa suitability: opportunity-cost (Zabel et al., 2024) ----
# cost.tif is cocoa crop suitability, clipped to the UGF at its native ~963m
# and rescaled to 0-1, where 1 is the most cocoa-suitable land in the region.
cost_path <- here("data/cocoa suitability/cost.tif")
if (!file.exists(cost_path)) {
  cocoa_native <- rast(here("data/cocoa suitability", "cocoa_suitability_ugf.tif"))
  # N/A cells are treated as 0 (no cocoa suitability)
  envelope <- rasterize(ugf_vect_3857, cocoa_native, field = 1, touches = TRUE)
  cocoa_native <- ifel(is.na(cocoa_native) & !is.na(envelope), 0, cocoa_native)
  cocoa_native <- cocoa_native / global(cocoa_native, "max", na.rm = TRUE)$max
  writeRaster(cocoa_native, cost_path, overwrite = TRUE)
  stamp("cost.tif absent -- rebuilt from cocoa_suitability_ugf.tif")
}

# ---- bring cocoa to the 3km grid ----
# Averaging, so a 3km cell carries the mean suitability of the ~963m cells
# inside it.
cost_native <- rast(cost_path)
cocoa_3km <- project(cost_native, template, method = "average")
# Same gap treatment as above, now at 3km.
boundary_rast <- rasterize(ugf_vect_3857, template, field = 1, touches = TRUE)
cocoa_norm <- ifel(is.na(cocoa_3km) & !is.na(boundary_rast), 0, cocoa_3km)
cocoa_norm <- mask(cocoa_norm, ugf_vect_3857, touches = TRUE)
cocoa_norm <- cocoa_norm / global(cocoa_norm, "max", na.rm = TRUE)$max

cost_uniform <- ifel(!is.na(cocoa_norm), 1, NA)
total_cells <- global(!is.na(cost_uniform), "sum", na.rm = TRUE)$sum
budget <- round(total_cells * 0.3)
stamp(sprintf("grid %.1f x %.1f m | planning units %d | budget (30%%) %d",
      res(template)[1], res(template)[2], total_cells, budget))

# ---- established protected areas ----
# Locked into every solution.
protected_areas <- st_read(here("data/protected areas", "protected_areas_ugf.shp"), quiet = TRUE)
protected_areas_rast <- rasterize(vect(protected_areas), cocoa_norm, field = 1, background = 0)
protected_areas_rast <- mask(protected_areas_rast, ugf_vect_3857, touches = TRUE)
stamp(sprintf("protected areas: %d cells (%.1f%% of PUs)",
      global(protected_areas_rast, "sum", na.rm = TRUE)$sum,
      100 * global(protected_areas_rast, "sum", na.rm = TRUE)$sum / total_cells))

# ---- urban land ----
# Locked out of the solution, except where it already sits inside a protected
# area -- those cells are locked IN, and a cell cannot be both.
urban <- rast(here("data/LC/urban_ugf.tif"))
urban_3km <- mask(project(urban, cocoa_norm, method = "near"), ugf_vect_3857, touches = TRUE)
urban_not_protected <- ifel(protected_areas_rast == 1, NA, urban_3km)

# ---- human pressure (gHM) ----
ghm_3km <- mask(project(rast(here("data/gHM/gHM.tif")), cocoa_norm, method = "average"),
                ugf_vect_3857, touches = TRUE)

# ---- species: 18 birds, Area of Habitat ----
# One layer per species: 1 where the species has habitat, NA elsewhere.
# Building the stack is slow, so it is cached to disk and reused -- but only
# after checking the cache still has the right number of layers and the right
# geometry, since a stale cache would silently poison every result below.
birds_star_t <- read.csv(here("data/birds_species_traits/birds_star_t.csv"))
clipped_dir <- here("data/birds_aoh/clipped")
aoh_cache <- here("data/birds_aoh/aoh_stack_3km.tif")

rebuild <- TRUE
if (file.exists(aoh_cache)) {
  cached <- rast(aoh_cache)
  if (nlyr(cached) == nrow(birds_star_t) && compareGeom(cached[[1]], cocoa_norm, stopOnError = FALSE)) {
    species_stack <- cached
    names(species_stack) <- birds_star_t$scientific_name
    rebuild <- FALSE
    stamp("loaded cached AOH stack (geometry verified)")
  }
}
if (rebuild) {
  species_stack <- NULL
  for (sp in birds_star_t$scientific_name) {
    f <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), "_R.tif"))
    if (!file.exists(f)) f <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), ".tif"))
    stopifnot(file.exists(f))
    r <- project(rast(f), cocoa_norm, method = "near")
    # These files store no-data as NaN inside an unsigned integer type, which
    # cannot represent NaN -- reprojecting turns those cells into 2^32 rather
    # than NA. Real habitat is always exactly 1, so treat anything else as
    # no-data.
    r <- mask(ifel(r == 1, 1, NA), ugf_vect_3857, touches = TRUE)
    names(r) <- sp
    species_stack <- if (is.null(species_stack)) r else c(species_stack, r)
    stamp(paste("  built", sp))
  }
  writeRaster(species_stack, aoh_cache, overwrite = TRUE)
}
species_stack_norm <- species_stack / global(species_stack, "max", na.rm = TRUE)$max

# ---- guard against missing gHM ----
# prioritizr reads an NA feature value as zero, so a planning unit with missing
# gHM would silently lose all of its habitat value. Report how many cells are
# affected, then treat missing gHM inside the study area as 0 (no measured
# human modification).
any_aoh <- !is.na(sum(species_stack_norm, na.rm = TRUE))
stamp(sprintf("PUs with NA gHM: %d | of which hold AOH: %d",
      global(!is.na(cost_uniform) & is.na(ghm_3km), "sum", na.rm = TRUE)$sum,
      global(any_aoh & is.na(ghm_3km), "sum", na.rm = TRUE)$sum))

ghm_pen   <- ifel(is.na(ghm_3km) & !is.na(cost_uniform), 0, ghm_3km)
cocoa_pen <- cocoa_norm

# ---- apply the habitat-quality discount ----
# Each species layer is scaled by (1 - gHM): habitat in a pristine cell counts
# fully, habitat in a heavily modified cell counts for little.
species_stack_eff <- species_stack_norm * (1 - ghm_pen)
names(species_stack_eff) <- names(species_stack_norm)
stamp(sprintf("raw AOH cell-sum %.0f -> discounted %.0f",
      sum(global(species_stack_norm, "sum", na.rm = TRUE)$sum),
      sum(global(species_stack_eff, "sum", na.rm = TRUE)$sum)))

# ---- save the layers ----
# Write the 3km layers out so they can be looked at, or reused, without
# re-running the whole script.
LAYERS <- here("outputs/layers_3km")
dir.create(LAYERS, recursive = TRUE, showWarnings = FALSE)
writeRaster(cocoa_norm,           file.path(LAYERS, "cocoa.tif"),             overwrite = TRUE)
writeRaster(ghm_pen,              file.path(LAYERS, "ghm.tif"),               overwrite = TRUE)
writeRaster(protected_areas_rast, file.path(LAYERS, "protected_areas.tif"),   overwrite = TRUE)
writeRaster(urban_not_protected,  file.path(LAYERS, "urban_locked_out.tif"),  overwrite = TRUE)
writeRaster(cost_uniform,         file.path(LAYERS, "planning_units.tif"),    overwrite = TRUE)
writeRaster(species_stack_eff,    file.path(LAYERS, "species_discounted.tif"), overwrite = TRUE)
stamp(sprintf("layers written to %s", LAYERS))

# ---- targets and weights ----
# Every species needs half its habitat covered. Shortfalls are weighted by IUCN
# category, so missing the target for a critically endangered species costs 64x
# what missing it for a least-concern species does. These are the same weights
# the CAPTAIN model uses, so the two approaches are penalising like for like.
species_targets_2 <- 0.5
species_traits <- read.csv(here("data/birds_species_traits/species_traits.csv"))
captain_iucn_weight <- c("1" = 1, "2" = 8, "3" = 16, "4" = 32, "5" = 64)
weights_iucn <- captain_iucn_weight[as.character(
  species_traits$conservation_status[match(names(species_stack), species_traits$species)])]
stopifnot(!anyNA(weights_iucn))


# ---- calibrating the cocoa penalty ----
# The penalty is set as a share of the shortfall you would get without the cocoa
# penalty. This step builds the problem without the cocoa penalty, solves it, 
# and records its shortfall and how much suitable cocoa land it selects.
p0 <- problem(cost_uniform, features = species_stack_eff) |>
  add_min_shortfall_objective(budget) |>
  add_relative_targets(species_targets_2) |>
  add_feature_weights(weights_iucn) |>
  add_locked_out_constraints(urban_not_protected) |>
  add_locked_in_constraints(protected_areas_rast) |>
  add_binary_decisions() |>
  add_gurobi_solver(gap = 0.001, verbose = FALSE)

weighted_shortfall <- function(p, s) sum(weights_iucn * eval_target_coverage_summary(p, s)$relative_shortfall)

s0 <- solve(p0)
obj0   <- weighted_shortfall(p0, s0)
cocoa0 <- global(cocoa_pen * s0, "sum", na.rm = TRUE)$sum
stamp(sprintf("calibration: shortfall %.3f | cells %d | sum(cocoa) %.1f",
      obj0, global(s0, "sum", na.rm = TRUE)$sum, cocoa0))

# Share of that baseline shortfall worth giving up to avoid cocoa land. Past
# 0.5, each further unit of coverage given up buys much less cocoa avoided.
COCOA_SHARE   <- 0.5
COCOA_PENALTY <- COCOA_SHARE * obj0 / cocoa0
stamp(sprintf("cocoa penalty: share %.2f -> %.6g", COCOA_SHARE, COCOA_PENALTY))

# ---- UGF prioritization problem ----
# The full problem, identical to the calibration problem above except for the
# add_linear_penalties() line.
p <- problem(cost_uniform, features = species_stack_eff) |>
  add_min_shortfall_objective(budget) |>
  add_relative_targets(species_targets_2) |>
  add_feature_weights(weights_iucn) |>
  add_linear_penalties(penalty = COCOA_PENALTY, data = cocoa_pen) |>
  add_locked_out_constraints(urban_not_protected) |>
  add_locked_in_constraints(protected_areas_rast) |>
  add_binary_decisions() |>
  add_gurobi_solver(gap = 0.001, verbose = FALSE)

s_gd <- solve(p)
writeRaster(s_gd, here("outputs/prioritizr_p3a_ghm_discount_3km.tif"), overwrite = TRUE)

n_sel <- global(s_gd, "sum", na.rm = TRUE)$sum
stamp(sprintf("solved: %d cells (%.1f%% of area) | weighted shortfall %.3f",
      n_sel, 100 * n_sel / total_cells, weighted_shortfall(p, s_gd)))
cat(sprintf("\nmean cocoa in selected cells %.4f (region-wide %.4f)\n",
    global(cocoa_pen * s_gd, "sum", na.rm = TRUE)$sum / n_sel,
    global(cocoa_pen, "mean", na.rm = TRUE)$mean))
cat(sprintf("mean gHM   in selected cells %.4f (region-wide %.4f)\n",
    global(ghm_pen * s_gd, "sum", na.rm = TRUE)$sum / n_sel,
    global(ghm_pen, "mean", na.rm = TRUE)$mean))

# ---- how much of each species' habitat was captured? ----
# 1) Raw: counts every habitat cell equally.
# 2) Quality-weighted: counts each cell by (1 - gHM).
q <- 1 - ghm_pen
coverage <- do.call(rbind, lapply(1:nlyr(species_stack_norm), function(i) {
  a <- !is.na(species_stack_norm[[i]])
  tot_raw  <- global(a, "sum", na.rm = TRUE)$sum
  tot_qual <- global(a * q, "sum", na.rm = TRUE)$sum
  f <- function(wt, tot) 100 * global(a * wt * (s_gd == 1), "sum", na.rm = TRUE)$sum / tot
  data.frame(species = names(species_stack_norm)[i], iucn_weight = weights_iucn[i],
             aoh_cells = tot_raw, pct_raw = round(f(1, tot_raw), 1),
             pct_qual = round(f(q, tot_qual), 1))
}))
cat("\n=== AOH coverage, 3km ===\n")
print(coverage[order(-coverage$iucn_weight, coverage$species), ], row.names = FALSE)
write.csv(coverage, here("outputs/p3a_ghm_discount_3km_coverage.csv"), row.names = FALSE)

big <- coverage[coverage$aoh_cells > 1000, ]
cat("\nIUCN-weighted mean % AOH captured\n")
cat(sprintf("  all 18 spp       raw %.1f | quality-weighted %.1f\n",
    weighted.mean(coverage$pct_raw, coverage$iucn_weight),
    weighted.mean(coverage$pct_qual, coverage$iucn_weight)))
cat(sprintf("  %2d spp >1000 px  raw %.1f | quality-weighted %.1f\n", nrow(big),
    weighted.mean(big$pct_raw, big$iucn_weight),
    weighted.mean(big$pct_qual, big$iucn_weight)))
stamp("solve stage complete")

############################################################################
# FIGURE
############################################################################
# Map of the solution. Colors match the CAPTAIN figures so the two sets can be
# compared side by side:
#   grey = not selected, blue = newly selected, dark grey = protected area

RES <- 150
LINE_PX <- 0.2 * RES        # one margin "line" in px at this resolution
MAR_SIDE <- 3.1
MAR_TOP <- 8.4
TITLE_LINE <- 5.6           # gap between map and title, in margin lines

fig_geom <- function(r, width = 2200, legend_px = 95) {
  e <- ext(r)
  asp <- (e$xmax - e$xmin) / (e$ymax - e$ymin)
  plot_w <- width - 2 * MAR_SIDE * LINE_PX
  panel_h <- plot_w / asp + MAR_TOP * LINE_PX
  list(width = width, height = round(panel_h + legend_px), heights = c(panel_h, legend_px))
}

# 0 not selected, 1 newly selected, 2 already protected
sel_class <- mask(ifel(protected_areas_rast == 1, 2, ifel(s_gd == 1, 1, 0)), s_gd)
sol_cols <- c("grey92", "#1f78b4", "grey40")

g <- fig_geom(sel_class, legend_px = 95)
png(here("outputs/solution_3a_ghm_discount_3km.png"),
    width = g$width, height = g$height, res = RES)
layout(matrix(c(1, 2), nrow = 2), heights = g$heights)
par(mar = c(0, MAR_SIDE, MAR_TOP, MAR_SIDE))
plot(sel_class, col = sol_cols, breaks = c(-0.5, 0.5, 1.5, 2.5), main = "", legend = FALSE)
mtext("PrioritizR solution (gHM-discounted AOH + cocoa penalty, 3km)",
      side = 3, line = TITLE_LINE, cex = 2.0, font = 2)
par(mar = c(0, 0, 0, 0))
plot.new()
legend(x = 0.5, y = 1, xjust = 0.5, yjust = 1,
       legend = c("Not selected", "Newly selected", "Protected Areas (locked-in)"),
       fill = sol_cols, ncol = 3, bty = "n", cex = 1.6, xpd = TRUE)
dev.off()

fr <- freq(sel_class); c0 <- function(v) { r <- fr$count[fr$value == v]; if (length(r)) r else 0 }
cat(sprintf("[3km] not selected %d | newly selected %d | protected areas %d\n",
    c0(0), c0(1), c0(2)))
stamp("figure written")
