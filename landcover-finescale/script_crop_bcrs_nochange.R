## Crop BCRs out of NLCD land cover data - NOT changeproduct

#### Libraries ####
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(stringr)

#### Data needed for all time steps ####
bcrs <- c(9, 12, 13, 14, 18, 19, 23, 27, 29) # BCRs of interest

# BCR shapefile
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us.proj <- readOGR("BCRs_contiguous_us.shp")

#### 1992 ####

## Crop raster data for each BCR
#setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_1992_landcover_2018_08_31/")
#filepath <- list.files(pattern = "tif$")
#region <- raster(filepath)

#crs.nlcd <- crs(region)
#bcrs.proj <- spTransform(us.proj, crs.nlcd)

# output: raster file of landcover change
# output: dataframe of landcover change of interest
#setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_92/")
#for(i in 1:length(bcrs)) {
#  bcr <- bcrs[i]
#  bcr.sub <- us.proj[us.proj@data$BCR == bcr, ]
#  bcr.sub1 <- bcr.sub[, -(7:8)]
#  area.sub <- crop(region, extent(bcr.sub1))
#  zones.sub <- raster::mask(area.sub, bcr.sub1)
#  filename <- paste0("nlcd_30x30_1992_bcr_", bcr, ".grd")
#  writeRaster(zones.sub, filename = filename)
#}

#### 2001 ####

## Crop raster data for each BCR
#setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_2001_landcover_2011_edition_2014_10_10/")
#filepath <- list.files(pattern = "img$")
#region <- raster(filepath)

#crs.nlcd <- crs(region)
#bcrs.proj <- spTransform(us.proj, crs.nlcd)

# output: raster file of landcover change
# output: dataframe of landcover change of interest
#setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_01/")
#for(i in 1:length(bcrs)) {
#  bcr <- bcrs[i]
#  bcr.sub <- us.proj[us.proj@data$BCR == bcr, ]
#  bcr.sub1 <- bcr.sub[, -(7:8)]
#  area.sub <- crop(region, extent(bcr.sub1))
#  zones.sub <- raster::mask(area.sub, bcr.sub1)
#  filename <- paste0("nlcd_30x30_2001_bcr_", bcr, ".grd")
#  writeRaster(zones.sub, filename = filename)
#}

#### 2006 ####

## Crop raster data for each BCR
setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_2006_landcover_2011_edition_2014_10_10/")
filepath <- list.files(pattern = "img$")
region <- raster(filepath)

crs.nlcd <- crs(region)
bcrs.proj <- spTransform(us.proj, crs.nlcd)

# output: raster file of landcover change
# output: dataframe of landcover change of interest
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_06/")
for(i in 2:length(bcrs)) {
  bcr <- bcrs[i]
  bcr.sub <- us.proj[us.proj@data$BCR == bcr, ]
  bcr.sub1 <- bcr.sub[, -(7:8)]
  area.sub <- crop(region, extent(bcr.sub1))
  zones.sub <- raster::mask(area.sub, bcr.sub1)
  filename <- paste0("nlcd_30x30_2006_bcr_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
}

#### 2011 ####

## Crop raster data for each BCR
setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_2011_landcover_2011_edition_2014_10_10/")
filepath <- list.files(pattern = "img$")
region <- raster(filepath)

crs.nlcd <- crs(region)
bcrs.proj <- spTransform(us.proj, crs.nlcd)

# output: raster file of landcover change
# output: dataframe of landcover change of interest
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_11/")
for(i in 1:length(bcrs)) {
  bcr <- bcrs[i]
  bcr.sub <- us.proj[us.proj@data$BCR == bcr, ]
  bcr.sub1 <- bcr.sub[, -(7:8)]
  area.sub <- crop(region, extent(bcr.sub1))
  zones.sub <- raster::mask(area.sub, bcr.sub1)
  filename <- paste0("nlcd_30x30_2011_bcr_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
}

