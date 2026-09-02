# Summary tables for the three PrioritizR-vs-CAPTAIN comparisons, using the
# same inputs and column definitions as run_characteristics.R.
#
#   Rscript scripts/prioritzr_captain3_comparison_table.R
#
# Table 1: agreement, on NEWLY selected cells only. Locked-in PAs are excluded
#   throughout -- both models contain the same 8,451 PA cells, which neither
#   chose, so counting them inflates agreement (per captain3_priority_to_raster.R).
#
# Table 2: what each model picked that the other did not -- mean cost,
#   disturbance and per-cell richness of the disagreement cells. This is where
#   the two models actually differ, so it is more diagnostic than comparing
#   their full selections.
#
# CAVEAT: the two models do not share a planning universe. CAPTAIN's 30% is of
# aoh_union_mask (44,004 cells -> 13,201); prioritizr's is of every costed cell
# inside the UGF (46,172 -> 13,852). Prioritizr therefore has 651 more cells and
# can select AOH-free land CAPTAIN could never reach. Differences below include
# that asymmetry.
suppressMessages({library(here); library(terra)})

DD  <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected"
ED  <- file.path(DD, "environmental_layers")
CLS <- c("LC", "NT", "VU", "EN", "CR")

tr   <- read.csv(file.path(DD, "species_traits.csv"))
cost <- as.matrix(rast(file.path(ED, "costs.tif")), wide = TRUE)
dist <- as.matrix(rast(file.path(ED, "disturbance.tif")), wide = TRUE)
dfut <- as.matrix(rast(file.path(ED, "disturbance_future_2075.tif")), wide = TRUE)
pa   <- as.matrix(rast(file.path(ED, "protected_areas.tif")), wide = TRUE); pa[is.na(pa)] <- 0
aoh  <- as.matrix(rast(here("outputs/prioritizr_p3a_ghm_discount_3km.tif")), wide = TRUE)

occ <- setNames(vector("list", length(CLS)), CLS)
for (i in seq_len(nrow(tr))) {
  s <- as.matrix(rast(file.path(DD, "present_habitat_suitability",
                                paste0(tr$species[i], ".tif"))), wide = TRUE) >= 0.5
  s[is.na(s)] <- FALSE
  k <- CLS[tr$conservation_status[i]]
  occ[[k]] <- if (is.null(occ[[k]])) s * 1 else occ[[k]] + s * 1
}
keep <- CLS[!vapply(occ, is.null, TRUE)]   # LC is empty: no LC species in the pool

pz <- as.matrix(rast(here("outputs/prioritizr_p3a_ghm_discount_3km.tif")), wide = TRUE)
pz[is.na(pz)] <- 0
pz_new <- pz > 0 & pa != 1

variants <- c(
  "CAPTAIN (simple)"                   = "outputs/captain3_seed232373165_corrected",
  "CAPTAIN (future)"                   = "outputs/captain3_seed232373165_corrected_future",
  "CAPTAIN (future + multi time step)" = "outputs/captain3_seed232373165_corrected_future_staged")

chars <- function(sel) c(cells = sum(sel),
  cost = mean(cost[sel], na.rm = TRUE),
  disturbance = mean(dist[sel], na.rm = TRUE),
  dist_2075 = mean(dfut[sel], na.rm = TRUE),
  vapply(keep, function(k) mean(occ[[k]][sel], na.rm = TRUE), 0))

agree_rows <- list(); diff_rows <- list()
for (nm in names(variants)) {
  g <- as.matrix(read.csv(here(variants[[nm]], "selected_grid.csv"), header = FALSE))
  cap_new <- g > 0 & pa != 1

  both   <- pz_new & cap_new
  only_p <- pz_new & !cap_new
  only_c <- cap_new & !pz_new

  agree_rows[[nm]] <- c(
    prioritizr = sum(pz_new), captain = sum(cap_new),
    both = sum(both), prioritizr_only = sum(only_p), captain_only = sum(only_c),
    IoU = sum(both) / sum(pz_new | cap_new),
    pct_captain_matched = 100 * sum(both) / sum(cap_new))

  diff_rows[[paste0(nm, "  >> PrioritizR only")]] <- chars(only_p)
  diff_rows[[paste0(nm, "  >> CAPTAIN only")]]    <- chars(only_c)
  diff_rows[[paste0(nm, "  >> both")]]            <- chars(both)
}

t1 <- do.call(rbind, agree_rows)
t2 <- do.call(rbind, diff_rows)

cat("\n=== Table 1: agreement on newly selected cells (PAs excluded) ===\n")
print(round(t1, 3))
cat("\n=== Table 2: characteristics of agreement / disagreement cells ===\n")
cat("species pool:", paste(sprintf("%s=%d", keep,
    vapply(keep, function(k) sum(tr$conservation_status == which(CLS == k)), 0L)), collapse = "  "), "\n\n")
print(round(t2, 3))

write.csv(t1, here("outputs/prioritzr_captain3_agreement_table.csv"))
write.csv(t2, here("outputs/prioritzr_captain3_disagreement_characteristics.csv"))
cat("\nwritten: outputs/prioritzr_captain3_agreement_table.csv\n")
cat("written: outputs/prioritzr_captain3_disagreement_characteristics.csv\n")
