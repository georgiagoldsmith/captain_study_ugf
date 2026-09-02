# Agreement: PrioritizR vs. CAPTAIN
# --------------------------------------

#   Rscript scripts/agreement_maps_captain_prioritizr.R

# Compares the PrioritizR solution against CAPTAIN, for each of the three
# CAPTAIN variants (all seed 232373165):
#
#   simple          present-day data, one time step
#   future          future data, one time step
#   future_staged   future data, protection spread over several time steps
#
# Produces one map per variant, plus a table of how much they overlap.
# Cells are classed as: neither, PrioritizR only, CAPTAIN only, both, or
# already protected.
#
# Protected areas are excluded from the overlap score. Both models contain all
# of them because they are locked in, so counting them would make the two look
# far more alike than they are.
#
# CAVEAT: the two models are not choosing from the same pool of land. CAPTAIN
# picks 30% of the cells that hold bird habitat (44,004 cells); PrioritizR picks
# 30% of every costed cell in the UGF (46,172). PrioritizR therefore has 651
# more cells to choose from, and can pick land with no bird habitat that CAPTAIN
# could never reach. Some of the disagreement below is that, not a difference
# of opinion.

suppressMessages({library(here); library(terra); library(sf)})

prz <- rast(here("outputs/prioritizr_p3a_ghm_discount_3km.tif"))
ugf <- vect(st_transform(st_read(here("data/UGF_gp.shp", "UGF_gp.shp"), quiet = TRUE), 3857))
pas <- rast(here("outputs/layers_3km/protected_areas.tif"))
stopifnot(compareGeom(pas, prz, stopOnError = FALSE))
pas <- mask(pas, ugf, touches = TRUE)

# CAPTAIN writes its plans as a plain matrix of 0/1, not a raster, so they are
# laid back onto the PrioritizR grid here.
read_captain <- function(dir) {
  m <- as.matrix(read.csv(file.path(here("outputs", dir), "selected_grid.csv"), header = FALSE))
  stopifnot(nrow(m) == nrow(prz), ncol(m) == ncol(prz))
  r <- rast(m); ext(r) <- ext(prz); crs(r) <- crs(prz)
  mask(r, prz)
}

variants <- list(
  list(dir = "captain3_seed232373165_corrected",               par = "(simple)",                   tag = "simple"),
  list(dir = "captain3_seed232373165_corrected_future",        par = "(future)",                   tag = "future"),
  list(dir = "captain3_seed232373165_corrected_future_staged", par = "(future + multi time step)", tag = "future_staged")
)

cols <- c("grey93", "#e66101", "#2c7bb6", "#1a9641", "grey40")
fmt <- function(n) format(n, big.mark = ",")

rows <- list()
for (v in variants) {
  cap <- read_captain(v$dir)
  cat(sprintf("\n=== %s ===\n", v$tag))

  pa <- pas == 1
  p_new <- (prz == 1) & !pa
  c_new <- (cap > 0.5) & !pa
  cls <- mask(ifel(pa, 4, ifel(p_new & c_new, 3, ifel(c_new, 2, ifel(p_new, 1, 0)))), prz)

  n <- freq(cls); cnt <- function(x) { r <- n$count[n$value == x]; if (length(r)) r else 0 }
  both <- cnt(3); only_p <- cnt(1); only_c <- cnt(2)
  iou <- both / (both + only_p + only_c)
  cat(sprintf("  both %d | prioritizr only %d | CAPTAIN only %d | PAs %d | IoU %.4f\n",
      both, only_p, only_c, cnt(4), iou))

  rows[[v$tag]] <- c(prioritizr = only_p + both, captain = only_c + both,
                     both = both, prioritizr_only = only_p, captain_only = only_c,
                     IoU = iou, pct_captain_matched = 100 * both / (only_c + both))

  # legend() sizes every column to the widest label, so the long
  # "(future + multi time step)" entry would overflow a 4-column row. Drop to
  # two columns when the label gets long.
  labs <- as.expression(list(
    bquote("neither"),
    bquote(.(sprintf("both (%s)", fmt(both)))),
    bquote(.(sprintf("PrioritizR (%s)", fmt(only_p)))),
    bquote("CAPTAIN" ~ italic(.(v$par)) ~ .(sprintf("(%s)", fmt(only_c))))))
  leg_cols <- c("grey93", "#1a9641", "#e66101", "#2c7bb6")
  ncol_leg <- if (nchar(v$par) > 12) 2 else 4
  leg_h <- if (ncol_leg == 2) 1.6 else 1.15

  png(here(sprintf("outputs/prioritzr_vs_captain3_%s.png", v$tag)),
      width = 2200, height = if (ncol_leg == 2) 1320 else 1250, res = 150)
  layout(matrix(1:2, nrow = 2), heights = c(6, leg_h))
  par(mar = c(0.2, 2.2, 6.5, 2.2))
  plot(cls, col = cols, breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5),
       legend = FALSE, main = "", axes = TRUE)
  title(main = bquote("PrioritizR vs CAPTAIN" ~ italic(.(v$par)) ~
                      .(sprintf("(IoU = %.2f)", iou))),
        cex.main = 1.85, line = 4.2)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend(x = 0.5, y = 0.66, xjust = 0.5, yjust = 0.5, legend = labs,
         fill = leg_cols, ncol = ncol_leg, bty = "n", cex = 1.3, xpd = TRUE)
  legend(x = 0.5, y = 0.18, xjust = 0.5, yjust = 0.5,
         legend = "Protected Areas (locked in)", fill = "grey40",
         ncol = 1, bty = "n", cex = 1.3, xpd = TRUE)
  dev.off()
  writeRaster(cls, here(sprintf("outputs/prioritzr_vs_captain3_%s.tif", v$tag)), overwrite = TRUE)
}

t1 <- do.call(rbind, rows)
cat("\n=== overlap on newly selected cells (protected areas excluded) ===\n")
print(round(t1, 3))
write.csv(t1, here("outputs/prioritzr_captain3_agreement_table.csv"))
cat("\nwritten: outputs/prioritzr_captain3_agreement_table.csv\n")
