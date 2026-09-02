# Characteristics of the cells a solution selected: mean cost, mean present and
# 2075 disturbance, and mean species richness per cell split by IUCN class.
#
#   Rscript scripts/run_characteristics.R                     # standard table:
#                                                             # 3 CAPTAIN variants
#                                                             # + prioritizr
#   Rscript scripts/run_characteristics.R <run_dir> [<run_dir> ...]
#
# Merged from run_characteristics.R and run_characteristics_with_prioritizr.R.
# Those two shared their entire setup, inputs and column definitions and
# differed only in which rows they built, so the difference is now the argument
# list rather than a second file. Both output filenames are unchanged.
#
# Richness is PRESENT-DAY occupancy (SDM >= MIN_SUITABILITY 0.5), i.e. how many
# species of each class actually live in an average selected cell -- not the
# fractional end-of-episode class totals in summary.txt, which are averages over
# the 30 stochastic replicates and so carry decimals for a different reason.
# Stats cover NEWLY selected cells only; the locked-in PAs are shared by every
# run and would wash out the differences between them.
#
# reward is CAPTAIN-internal (final-epoch avg_reward from its RL training log).
# prioritizr has no episode and no reward, so that cell is NA by construction,
# not missing data.
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
  outfile <- "outputs/run_characteristics.csv"
} else {
  # ---- standard table: the three CAPTAIN variants, then the ILP ----
  captain_runs <- c(
    "CAPTAIN (simple)"                   = "outputs/captain3_seed232373165_corrected",
    "CAPTAIN (future)"                   = "outputs/captain3_seed232373165_corrected_future",
    "CAPTAIN (future + multi time step)" = "outputs/captain3_seed232373165_corrected_future_staged")
  for (nm in names(captain_runs)) rows[[nm]] <- captain_row(captain_runs[[nm]])

  pz <- as.matrix(rast(here("outputs/prioritizr_p3a_ghm_discount_3km.tif")), wide = TRUE)
  pz[is.na(pz)] <- 0
  rows[["PrioritizR (gHM discount + cocoa penalty)"]] <- metrics(pz > 0 & pa != 1, NA_real_)

  outfile <- "outputs/run_characteristics_with_prioritizr.csv"
}

out <- do.call(rbind, rows)
cat("\nspecies pool by class:", paste(sprintf("%s=%d", CLS,
    vapply(CLS, function(k) sum(tr$conservation_status == which(CLS == k)), 0L)), collapse = "  "), "\n\n")
print(round(out, 3))
write.csv(out, here(outfile))
cat(sprintf("\nwritten: %s\n", outfile))
