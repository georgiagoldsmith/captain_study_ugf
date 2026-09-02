# Context map: UGF boundary shown against African country outlines.
# Country boundaries come from rnaturalearthdata (bundled locally with the
# rnaturalearth package -- no network download needed).

library(here)
library(sf)
library(rnaturalearth)

africa <- ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")

ugf_boundary <- st_read(here("data/UGF_gp.shp", "UGF_gp.shp"), quiet = TRUE)

# match CRS -- UGF boundary is already in Web Mercator (3857), reproject
# Africa to match rather than the other way around (keeps UGF's own
# geometry untouched)
africa_3857 <- st_transform(africa, st_crs(ugf_boundary))

# semi-transparent UGF fill (solid border keeps the boundary crisp)
ugf_fill <- adjustcolor("darkgreen", alpha.f = 0.5)

# size the map panel to the data's own aspect ratio (instead of guessing a
# heights ratio) so asp=1 doesn't leave blank space between the map and the
# legend panel below it
img_w <- 2400
africa_bbox <- st_bbox(africa_3857)
data_aspect <- (africa_bbox["xmax"] - africa_bbox["xmin"]) / (africa_bbox["ymax"] - africa_bbox["ymin"])
map_panel_h <- img_w / data_aspect
legend_panel_h <- 220

png(here("outputs/ugf_africa_context_map.png"), width = img_w,
    height = map_panel_h + legend_panel_h, res = 200)
layout(matrix(c(1, 2), nrow = 2), heights = c(map_panel_h, legend_panel_h))

par(mar = c(1.1, 1.1, 3.1, 1.1))
plot(st_geometry(africa_3857), col = "grey90", border = "grey50", lwd = 0.5,
     main = "Upper Guinean Forest (UGF) Study Area in West Africa")
plot(st_geometry(ugf_boundary), col = ugf_fill, border = "darkgreen", add = TRUE)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("African countries", "UGF study area"),
       fill = c("grey90", ugf_fill),
       border = c("grey50", "darkgreen"),
       ncol = 2, cex = 1.2, bty = "n", xpd = TRUE)

dev.off()

cat("Saved outputs/ugf_africa_context_map.png\n")

# ---- zoomed-in West Africa version ----
# buffer the UGF bbox by ~450km on each side so neighboring countries
# (Senegal to Nigeria/Cameroon) are visible for context
ugf_bbox <- st_bbox(ugf_boundary)
buffer_m <- 450000
# south of the coastline is open ocean with nothing to show -- a small
# buffer there is enough, vs. the full 450km used on the other 3 sides for
# neighboring-country context
buffer_south_m <- 60000
xlim <- c(ugf_bbox["xmin"] - buffer_m, ugf_bbox["xmax"] + buffer_m)
ylim <- c(ugf_bbox["ymin"] - buffer_south_m, ugf_bbox["ymax"] + buffer_m)

# country label points (centroids), clipped to the zoomed extent so labels
# for far-off countries don't get placed off-screen
africa_centroids <- st_centroid(africa_3857)
centroid_coords <- st_coordinates(africa_centroids)
in_view <- centroid_coords[, 1] > xlim[1] & centroid_coords[, 1] < xlim[2] &
           centroid_coords[, 2] > ylim[1] & centroid_coords[, 2] < ylim[2]

zoom_data_aspect <- diff(xlim) / diff(ylim)
zoom_map_panel_h <- img_w / zoom_data_aspect
zoom_legend_panel_h <- 220

png(here("outputs/ugf_africa_context_map_zoomed.png"), width = img_w,
    height = zoom_map_panel_h + zoom_legend_panel_h, res = 200)
layout(matrix(c(1, 2), nrow = 2), heights = c(zoom_map_panel_h, zoom_legend_panel_h))

par(mar = c(1.1, 1.1, 3.1, 1.1))
plot(st_geometry(africa_3857), col = "grey90", border = "grey50", lwd = 0.7,
     main = "Upper Guinean Forest (UGF) Study Area in West Africa",
     xlim = xlim, ylim = ylim)
plot(st_geometry(ugf_boundary), col = ugf_fill, border = "darkgreen", add = TRUE)

# white halo behind each label so names stay legible whether they land on
# the grey country background or the dark green UGF fill (e.g. Sierra
# Leone, Liberia)
label_x <- centroid_coords[in_view, 1]
label_y <- centroid_coords[in_view, 2]
label_txt <- africa$name[in_view]
halo_offset <- diff(xlim) * 0.0018
for (dx in c(-1, 0, 1)) {
  for (dy in c(-1, 0, 1)) {
    if (dx != 0 || dy != 0) {
      text(label_x + dx * halo_offset, label_y + dy * halo_offset,
           labels = label_txt, cex = 0.65, col = "white", font = 2)
    }
  }
}
text(label_x, label_y, labels = label_txt, cex = 0.65, col = "grey20", font = 2)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("African countries", "UGF study area"),
       fill = c("grey90", ugf_fill),
       border = c("grey50", "darkgreen"),
       ncol = 2, cex = 1.2, bty = "n", xpd = TRUE)

dev.off()

cat("Saved outputs/ugf_africa_context_map_zoomed.png\n")
