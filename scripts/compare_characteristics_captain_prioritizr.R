# PrioritizR & CAPTAIN solution characteristics
# --------------------------------------------------

#   Rscript scripts/compare_characteristics_captain_prioritizr.R               # standard table
#   Rscript scripts/compare_characteristics_captain_prioritizr.R <run_dir> ...  # runs you name

# One row per version, with the average cost, disturbance now, future 
# disturbance, and number of species per cell broken down by IUCN category. 
# Builds the comparison table: three CAPTAIN variants plus prioritizr.
#
# Species counts are how many species of each category actually live in an
# average chosen cell, counting a species as present where its habitat
# suitability is 0.5 or higher.
#
# Only newly chosen cells are counted. Protected areas are locked in, so 
# including them would inflate agreement.
#
# Reward comes from CAPTAIN's own training log. PrioritizR does not produce
# a reward, so that cell is empty.
#
# A second table describes which cells the models disagreed on:what PrioritizR
# picked that CAPTAIN did not, what CAPTAIN picked that PrioritizR did not, and
# what they both picked. 
suppressMessages({library(here); library(terra)})

DD  <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected"
ED  <- file.path(DD, "environmental_layers")
CLS <- c("LC", "NT", "VU", "EN", "CR")   # conservation_status 1..5

tr   <- read.csv(file.path(DD, "species_traits.csv"))
cost <- as.matrix(rast(file.path(ED, "costs.tif")), wide = TRUE)
dist <- as.matrix(rast(file.path(ED, "disturbance.tif")), wide = TRUE)
dfut <- as.matrix(rast(file.path(ED, "disturbance_future_2075.tif")), wide = TRUE)
pa   <- as.matrix(rast(file.path(ED, "protected_areas.tif")), wide = TRUE); pa[is.na(pa)] <- 0

# stack of per-species present-day occupancy, kept as one matrix per IUCN class
occ <- setNames(vector("list", length(CLS)), CLS)
for (i in seq_len(nrow(tr))) {
  s <- as.matrix(rast(file.path(DD, "present_habitat_suitability",
                                paste0(tr$species[i], ".tif"))), wide = TRUE) >= 0.5
  s[is.na(s)] <- FALSE
  k <- CLS[tr$conservation_status[i]]
  occ[[k]] <- if (is.null(occ[[k]])) s * 1 else occ[[k]] + s * 1
}

metrics <- function(sel, reward) {
  c(cells = sum(sel), reward = reward,
    cost = mean(cost[sel], na.rm = TRUE),
    disturbance = mean(dist[sel], na.rm = TRUE),
    dist_2075 = mean(dfut[sel], na.rm = TRUE),
    vapply(CLS, function(k) if (is.null(occ[[k]])) NA_real_ else mean(occ[[k]][sel], na.rm = TRUE), 0))
}

# one CAPTAIN run directory -> one row, reward pulled from its training log
captain_row <- function(rd) {
  g  <- as.matrix(read.csv(here(rd, "selected_grid.csv"), header = FALSE))
  lg <- readLines(here(rd, "training_log.tsv"))
  metrics(g > 0 & pa != 1, as.numeric(strsplit(lg[length(lg)], "\t")[[1]][3]))
}

args <- commandArgs(trailingOnly = TRUE)
rows <- list()

if (length(args)) {
  # ---- ad-hoc mode: one row per run directory named on the command line ----
  for (rd in args) {
    if (!file.exists(here(rd, "selected_grid.csv"))) { cat("missing:", rd, "\n"); next }
    rows[[rd]] <- captain_row(rd)
  }
  outfile <- "outputs/captain_run_characteristics.csv"
} else {
  # ---- standard table: the three CAPTAIN variants, then the ILP ----
  captain_runs <- c(
    "CAPTAIN (simple)"                   = "outputs/captain_seed232373165_simple",
    "CAPTAIN (future)"                   = "outputs/captain_seed232373165_future",
    "CAPTAIN (future + multi time step)" = "outputs/captain_seed232373165_future_staged")
  for (nm in names(captain_runs)) rows[[nm]] <- captain_row(captain_runs[[nm]])

  pz <- as.matrix(rast(here("outputs/prioritizr/prioritizr_solution.tif")), wide = TRUE)
  pz[is.na(pz)] <- 0
  rows[["PrioritizR (gHM discount + cocoa penalty)"]] <- metrics(pz > 0 & pa != 1, NA_real_)

  outfile <- "outputs/comparison/prioritizr_captain_characteristics.csv"
}

dir.create(here("outputs/comparison"), recursive = TRUE, showWarnings = FALSE)
out <- do.call(rbind, rows)
cat("\nspecies pool by class:", paste(sprintf("%s=%d", CLS,
    vapply(CLS, function(k) sum(tr$conservation_status == which(CLS == k)), 0L)), collapse = "  "), "\n\n")
print(round(out, 3))
write.csv(out, here(outfile))
cat(sprintf("\nwritten: %s\n", outfile))

# ---- second table: the cells the two models disagreed on ----
# Same columns as above, but each row is a slice of one comparison rather than a
# whole plan, so there is no reward to report.
if (!length(args)) {
  diff_rows <- list()
  for (nm in names(captain_runs)) {
    g <- as.matrix(read.csv(here(captain_runs[[nm]], "selected_grid.csv"), header = FALSE))
    cap_new <- g > 0 & pa != 1
    pz_new  <- pz > 0 & pa != 1
    diff_rows[[paste0(nm, "  >> PrioritizR only")]] <- metrics(pz_new & !cap_new, NA_real_)
    diff_rows[[paste0(nm, "  >> CAPTAIN only")]]    <- metrics(cap_new & !pz_new, NA_real_)
    diff_rows[[paste0(nm, "  >> both")]]            <- metrics(pz_new & cap_new,  NA_real_)
  }
  t2 <- do.call(rbind, diff_rows)
  t2 <- t2[, colnames(t2) != "reward"]
  cat("\n=== characteristics of the agreement / disagreement cells ===\n")
  print(round(t2, 3))
  write.csv(t2, here("outputs/comparison/prioritizr_captain_disagreement_characteristics.csv"))
  cat("\nwritten: outputs/comparison/prioritizr_captain_disagreement_characteristics.csv\n")
}
