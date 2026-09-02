library(sf)
library(here)
library(ggplot2)
library(terra)

##################################
# Protected Areas

## load data
protected_areas_0 <- st_read(here("data/protected areas/WDPA_WDOECM_Jun2026_Public_AF_shp/WDPA_WDOECM_Jun2026_Public_AF_shp_0", "WDPA_WDOECM_Jun2026_Public_AF_shp-polygons.shp"))

protected_areas_1 <- st_read(here("data/protected areas/WDPA_WDOECM_Jun2026_Public_AF_shp/WDPA_WDOECM_Jun2026_Public_AF_shp_1", "WDPA_WDOECM_Jun2026_Public_AF_shp-polygons.shp"))

protected_areas_2 <- st_read(here("data/protected areas/WDPA_WDOECM_Jun2026_Public_AF_shp/WDPA_WDOECM_Jun2026_Public_AF_shp_2", "WDPA_WDOECM_Jun2026_Public_AF_shp-polygons.shp"))

ugf_boundary <- st_read(here("data/UGF_gp.shp", "UGF_gp.shp"))

## join protected areas
protected_areas_ugf <- dplyr::bind_rows(
  protected_areas_0,
  protected_areas_1,
  protected_areas_2)

## clip
protected_areas_ugf <- st_transform(protected_areas_ugf, st_crs(ugf_boundary))
protected_areas_ugf <- st_intersection(protected_areas_ugf, ugf_boundary)

## remove UNESCO-MAB Biosphere Reserve

protected_areas_ugf <- st_read("data/protected areas/protected_areas_ugf.shp")
protected_areas_ugf <- protected_areas_ugf %>%
  filter(DESIG != "UNESCO-MAB Biosphere Reserve")

protected_areas_ugf <- st_intersection(protected_areas_ugf, ugf_boundary)

protected_areas_ugf <- st_make_valid(protected_areas_ugf)

protected_areas_ugf <- st_collection_extract(protected_areas_ugf, "POLYGON")

## plot -- styled to match the base-R terra maps (bold centered title sitting
## just above the plot box, white background, no gridlines, black border)
plot_pa_ugf <- ggplot() +
  geom_sf(data = protected_areas_ugf, aes(fill = DESIG), color = NA) +
  theme_bw() +
  labs(title = "Protected Areas in UGF", fill = "Designation") +
  theme(legend.position = "bottom") +
  theme(legend.key.size = unit(0.3, "cm"),  # shrink legend boxes
        legend.text = element_text(size = 10),  # shrink legend text
        legend.title = element_text(size = 12),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 24,
                                   margin = margin(b = 8)))
plot_pa_ugf

st_write(protected_areas_ugf, here("data/protected areas/protected_areas_ugf.shp"), append = FALSE)

(ggsave(here("outputs/protected_areas_ugf.png"), width = 17, height = 10, dpi = 300))
##########################################################################
#LC

library(terra)
library(sf)

# load original
lccs <- rast(here("data/LC", "C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.area-subset.15.10.-5.-18.nc"))

# transform boundary to match raster CRS (keeping raster in its native CRS)
ugf_native <- st_transform(ugf_boundary, crs(lccs))
ugf_vect_native <- vect(ugf_native)

# clip in native CRS
lccs_ugf_raw <- crop(lccs, ugf_vect_native)
lccs_ugf_raw <- mask(lccs_ugf_raw, ugf_vect_native)

# select layer
lccs_ugf <- lccs_ugf_raw["lccs_class"]

# convert to integer to preserve class values
lccs_ugf <- as.int(lccs_ugf)

# only reproject AFTER clipping, using nearest neighbor
lccs_ugf <- project(lccs_ugf, "EPSG:3857", method = "near")

# check
unique(lccs_ugf)

# save
writeRaster(lccs_ugf, here("data/LC/lccs_ugf.tif"), overwrite = TRUE)

lccs_ugf <- rast(here("data/LC/lccs_ugf.tif"))

# urban areas (class 190)
urban_ugf <- lccs_ugf == 190
writeRaster(urban_ugf, here("data/LC/urban_ugf.tif"), overwrite = TRUE)

# cropland (classes 10, 11, 20)
cropland_ugf <- ifel(is.na(lccs_ugf), NA, lccs_ugf %in% c(10, 11, 20))
writeRaster(cropland_ugf, here("data/LC/cropland_ugf.tif"), overwrite = TRUE)
plot(cropland_ugf)


###############################################################

# Cocoa crop suitability -> clipped to the UGF.
# This is the cocoa layer the analysis uses: cocoa_suitability_ugf.tif is
# normalised into cost.tif by prioritzr_ghm_discount_3km.R and read from there
# by the CAPTAIN export.
cocoa_suitability <- rast(here("data/cocoa suitability/cocoa_crop_suitability.tif"))
cocoa_suitability <- project(cocoa_suitability, crs(ugf_boundary))

cocoa_suitability_ugf <- crop(cocoa_suitability, ugf_boundary)
cocoa_suitability_ugf <- mask(cocoa_suitability_ugf, ugf_boundary)
plot(cocoa_suitability_ugf)
res(cocoa_suitability_ugf)

writeRaster(cocoa_suitability_ugf, here("data/cocoa suitability/cocoa_suitability_ugf.tif"), overwrite = TRUE)

###################################################################

#Roads

install.packages("osmdata")
library(osmdata)

# load and combine road shapefiles
roads_gh <- st_read(here("data/roads/ghana/gis_osm_roads_free_1.shp"))
roads_ci <- st_read(here("data/roads/ivory_coast/gis_osm_roads_free_1.shp"))
roads_gn <- st_read(here("data/roads/guinea/gis_osm_roads_free_1.shp"))
roads_sl <- st_read(here("data/roads/sierra_leone/gis_osm_roads_free_1.shp"))
roads_lr <- st_read(here("data/roads/liberia/gis_osm_roads_free_1.shp"))
roads_tg <- st_read(here("data/roads/togo/gis_osm_roads_free_1.shp"))
roads_bn <- st_read(here("data/roads/benin/gis_osm_roads_free_1.shp"))

# combine
roads_all <- bind_rows(roads_gh, roads_ci, roads_gn, roads_sl, roads_lr)

# filter to major roads only
roads_main <- roads_all |>
  filter(fclass %in% c("motorway", "trunk", "primary", "secondary"))

# clip to UGF
roads_main <- st_transform(roads_main, st_crs(ugf_boundary))
roads_ugf <- st_intersection(roads_main, ugf_boundary)


