# Crop BBS routes from BCRs

#### Libraries ####
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(stringr)

# Read in BBS routes shapefile

# 92-01
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/92-01/")
files <- list.files(pattern = "grd")
routes <- readOGR("/proj/hurlbertlab/bbs/BBS_routes_USandCanada/bbs_routes_usandcanada.shp")

for(i in 1:length(files)) {
  file <- files[i]
  bcr <- substr(file, 22, 27)
  bcr.raster <- raster(file)
  routes.sub <- crop(routes, extent(bcr.raster))
  zones.sub <- mask(bcr.raster, routes.sub)
  filename <- paste0("nlcd_30x30_1992_2001_routes_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
}

# 01-06
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/01-06/")
files <- list.files(pattern = "grd")

for(i in 1:length(files)) {
  file <- files[i]
  bcr <- substr(file, 22, 27)
  bcr.raster <- raster(file)
  area.sub <- crop(routes, extent(bcr.raster))
  zones.sub <- mask(area.sub, bcr.raster)
  filename <- paste0("nlcd_30x30_2001_2006_routes_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
}

# 06-11
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/06-11/")
files <- list.files(pattern = "grd")

for(i in 1:length(files)) {
  file <- files[i]
  bcr <- substr(file, 22, 27)
  bcr.raster <- raster(file)
  area.sub <- crop(routes, extent(bcr.raster))
  zones.sub <- mask(area.sub, bcr.raster)
  filename <- paste0("nlcd_30x30_2006_2011_routes_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
}
