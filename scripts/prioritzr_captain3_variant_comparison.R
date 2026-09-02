# PrioritizR (gHM-discounted AOH + cocoa penalty, 3km) vs CAPTAIN 3, across the
# three CAPTAIN model variants. Seed 232373165 is the only seed with all three
# runs present, so all three panels use it -- differences between panels are
# the model variant, not the seed.
#
#   _corrected                 simple: present-day data, one time step
#   _corrected_future          future data, one time step
#   _corrected_future_staged   future data, multiple time steps
#
# CAPTAIN runs are stored as selected_grid.csv matrices, not rasters, so they
# are read onto the 3km template the way captain3_priority_to_raster.R does.
# Locked-in PAs get their own class and are excluded from the IoU, per that
# script's convention: both models contain the same PA cells, which neither
# chose, so counting them as agreement inflates it.

library(here)
library(terra)
library(sf)

prz <- rast(here("outputs/prioritizr_p3a_ghm_discount_3km.tif"))
ugf <- vect(st_transform(st_read(here("data/UGF_gp.shp", "UGF_gp.shp"), quiet = TRUE), 3857))
pas <- rasterize(vect(st_read(here("data/protected areas", "protected_areas_ugf.shp"), quiet = TRUE)),
                 prz, field = 1, background = 0)
pas <- mask(pas, ugf, touches = TRUE)

read_captain <- function(dir) {
  m <- as.matrix(read.csv(file.path(here("outputs", dir), "selected_grid.csv"), header = FALSE))
  stopifnot(nrow(m) == nrow(prz), ncol(m) == ncol(prz))
  r <- rast(m); ext(r) <- ext(prz); crs(r) <- crs(prz)
  mask(r, prz)
}

variants <- list(
  list(dir = "captain3_seed232373165_corrected",               par = "(simple)",                      tag = "simple"),
  list(dir = "captain3_seed232373165_corrected_future",        par = "(future)",                      tag = "future"),
  list(dir = "captain3_seed232373165_corrected_future_staged", par = "(future + multi time step)",    tag = "future_staged")
)

# Styling matches plot_future_vs_static.R exactly so this reads as the same
# figure family as the CAPTAIN-vs-CAPTAIN comparisons:
#   grey93 neither / #e66101 first model / #2c7bb6 second model /
#   #1a9641 both / grey40 protected areas, lowercase labels with comma-grouped
#   counts, axes on, title via title(line = 4.2).
cols <- c("grey93", "#e66101", "#2c7bb6", "#1a9641", "grey40")
fmt <- function(n) format(n, big.mark = ",")

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

  # Legend order as requested: neither, both, PrioritizR, CAPTAIN (variant).
  # legend() sizes every column to the WIDEST label, so the long
  # "(future + multi time step)" entry would blow out a 4-column row -- fall
  # back to two columns when the labels get long rather than letting the box
  # overflow the device.
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
cat("\ndone\n")
