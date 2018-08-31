# Script: examine land cover change in fine scale NCLD data
# 1992-2006 for eastern deciduous and western grassland Bird Conservation Regions

#### Libraries ####
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(stringr)

#### Functions ####
# Get area changeproduct raster file from NLCD directory

# 2001-2011
#setwd("/proj/hurlbertlab/nlcd_landcover_change/")
#files <- list.files()
#nlcd.files <- files[str_detect(files, "2006")]
#dir.img <- getwd()

#get.file.img <- function(x) {
#  files2 <- list.files(paste0(dir.img, "/", nlcd.files[x], ""))
#  file.path <- files2[str_detect(files2, "img")]
#  return(list(folder = nlcd.files[x], file.name = file.path))
#}

#### Data needed for all time steps ####
bcrs <- c(9, 12, 13, 14, 18, 19, 23, 27, 29) # BCRs of interest

# BCR shapefile
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us.proj <- readOGR("BCRs_contiguous_us.shp")

#### 1992 - 2001 ####

## Crop raster data for each BCR
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/")
region <- raster("nlcd_1992_2001_changepixels_30m_30m.grd")

crs.nlcd <- crs(region)
bcrs.proj <- spTransform(us.proj, crs.nlcd)

# output: raster file of landcover change
# output: dataframe of landcover change of interest
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/")
for(i in 1:length(bcrs)) {
  bcr <- bcrs[i]
  bcr.sub <- us.proj[us.proj@data$BCR == bcr, ]
  area.sub <- crop(region, extent(bcr.sub))
  zones.sub <- extract(area.sub, bcr.sub, buffer = 200)
  filename <- paste0("nlcd_30x30_1992_2001_bcr_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
}

#### 2001 - 2011 ####

## 2001 - 2006

# Read in land cover data
#file.2001 <- get.file.img(1)
#nlcd2001 <- raster(paste0(dir.img, "/", file.2001$folder, "/", file.2001$file.name))

# BCRs
#us01 <- sp::spTransform(us.proj, crs(nlcd2001))

# output: raster file of landcover change
# output: dataframe of landcover change of interest
#setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/")

#for(i in 1:length(bcrs)) {
 # bcr <- bcrs[i]
 # bcr.sub <- us01[us01@data$BCR == bcr, ]
 # area.sub <- crop(nlcd2001, extent(bcr.sub))
 # zones.sub <- mask(area.sub, bcr.sub)
 # filename <- paste0("nlcd_30x30_2001_2006_bcr_", bcr, ".grd")
 # writeRaster(zones.sub, filename = filename, overwrite = T)
#}

## 2006-2011

# Read in land cover data
#file.2006 <- get.file.img(2)
#nlcd2006 <- raster(paste0(dir.img, "/", file.2006$folder, "/", file.2006$file.name))

# output: raster file of landcover change
# output: dataframe of landcover change of interest
#setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/")

#for(i in 1:length(bcrs)) {
#  bcr <- bcrs[i]
#  bcr.sub <- us01[us01@data$BCR == bcr, ]
#  area.sub <- crop(nlcd2006, extent(bcr.sub))
#  zones.sub <- mask(area.sub, bcr.sub)
#  filename <- paste0("nlcd_30x30_2006_2011_bcr_", bcr, ".grd")
#  writeRaster(zones.sub, filename = filename, overwrite = T)
#}
