# Calculating future (2075) disturbance in gHM units
# --------------------------------------------------

#   Rscript scripts/export_future_disturbance.R            # 3km corrected (default)
#   DATA_DIR=.../ugf_data_1km Rscript scripts/export_future_disturbance.R

# There is no future gHM product, so the SSP4-RCP3.4 LULC classes are calibrated
# against present-day gHM (median gHM of each 2020 class) and the class-to-class
# change is applied to the observed present gHM:
#   future = present_gHM + (median[class_2075] - median[class_2020])

# A cell whose land cover does not change gets zero delta.
# Urban and water in 2075 take the 1.0 exclusion, matching the present layer.
suppressMessages({library(here);library(terra);library(sf)})
DD <- Sys.getenv("DATA_DIR",
  "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected")
cat("target build:", basename(DD), "\n")
ugf  <- vect(st_transform(st_read(here("data/UGF_gp.shp","UGF_gp.shp"),quiet=TRUE),3857))
grid <- rast(file.path(DD,"environmental_layers","costs.tif"))
ghm  <- mask(project(rast(here("data/gHM/gHM.tif")), grid, method="average"), ugf, touches=TRUE)

get_lu <- function(year){
  r <- rast(here("data/future lulc/SSP4_RCP34", paste0("global_SSP4_RCP34_",year,".tif")))
  u <- project(ugf, crs(r))
  r <- mask(crop(r, ext(u)+0.05), u)
  mask(project(r, grid, method="near"), ugf, touches=TRUE)
}
lu20 <- get_lu(2020); lu75 <- get_lu(2075)

g <- values(ghm); l <- values(lu20); ok <- !is.na(g) & !is.na(l) & l>0
med <- sapply(1:7, function(k) if (sum(ok & l==k)>0) median(g[ok & l==k]) else NA)
cat("class medians (present gHM):", paste(sprintf("%d=%.3f",1:7,med), collapse=" "), "\n")

m20 <- classify(lu20, cbind(1:7, med)); m75 <- classify(lu75, cbind(1:7, med))
fut <- ghm + (m75 - m20)
fut <- clamp(fut, 0, 0.999)
fut <- ifel(lu75 %in% c(1,6), 1, fut)     # water, urban -> exclusion sentinel
writeRaster(fut, file.path(DD,"environmental_layers","disturbance_future_2075.tif"), overwrite=TRUE)
# keep a copy in the project too, next to the present-day gHM it is derived from
writeRaster(fut, here("data/gHM/disturbance_future_2075.tif"), overwrite=TRUE)

a <- values(ghm); b <- values(fut); k <- !is.na(a) & !is.na(b)
cat(sprintf("present mean %.3f -> future mean %.3f\n", mean(a[k]), mean(b[k])))
cat(sprintf("cells increasing %d | decreasing %d | unchanged %d\n",
    sum(b[k]>a[k]+1e-6), sum(b[k]<a[k]-1e-6), sum(abs(b[k]-a[k])<=1e-6)))
cat(sprintf("cells at the 1.0 sentinel: present %d -> future %d\n", sum(a[k]>=1), sum(b[k]>=1)))
