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

# 1992-2001
#setwd("nlcd_1992_to_2001_landcover_change")
#files <- list.files()
#area.files <- files[str_detect(files, "area_")]
#dir <- getwd()

#get.file <- function(x) {
#  files2 <- list.files(paste0(dir, "/", area.files[x], ""))
#  file.path <- files2[str_detect(files2, "area")]
#  return(list(folder = area.files[x], file.name = file.path))
#}

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

# Function to merge two raster files with different origins and extents, same resolution and crs
#merge.areas <- function(x, y) {
#  extent1 <- extent(x)
#  extent2 <- extent(y)
#  z <- raster(xmn = min(extent1[1], extent2[1]), xmx = max(extent1[2], extent2[2]), # make blank raster file with extent of the two files added
#              ymn = min(extent1[3], extent2[3]), ymx = max(extent1[4], extent2[4]), 
#              resolution = res(x), crs = crs(x))
#  merge1 <- mosaic(z, x, fun = max, tolerance = max(abs(origin(z) - origin(x)))) # function: where rasters overlap keep larger value
#  merge2 <- mosaic(merge1, y, fun = max, tolerance = max(abs(origin(merge1) - origin(y)))) # tolerance allows for different origins 
#  return(merge2) # returns the two rasters merged into one
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

#crs.nlcd <- crs(region)
#bcrs.proj <- spTransform(us.proj, crs.nlcd)

#setwd("/proj/hurlbertlab/gdicecco/NLCD_fragmentation/landcover-finescale/")
#codes <- read.csv("anderson_land_cover_codes.csv", stringsAsFactors = F) # NLCD land cover class codes

#fragcodes <- codes %>%
#  filter(grepl("Forest|Grassland", anderson))

# output: raster file of landcover change
# output: dataframe of landcover change of interest
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/")
for(i in 1:length(bcrs)) {
  bcr <- bcrs[i]
  bcr.sub <- bcrs.proj[bcrs.proj@data$BCR == bcr, ]
  area.sub <- crop(region, extent(bcr.sub))
  zones.sub <- mask(area.sub, bcr.sub)
  filename <- paste0("nlcd_30x30_1992_2001_bcr_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
#  area.df <- rasterToPoints(zones.sub, fun = fun(x) {x %in% fragcodes$modified})
#  write.csv(area.df, paste0("nlcd_30x30_1992_2001_bcr_", bcr, ".csv"), row.names = F)
}

#### 2001 - 2011 ####

## 2001 - 2006

# Read in land cover data
#file.2001 <- get.file.img(1)
#nlcd2001 <- raster(paste0(dir.img, "/", file.2001$folder, "/", file.2001$file.name))

# BCRs
#us01 <- sp::spTransform(us.proj, crs(nlcd2001))

# Land cover class codes
#setwd("/proj/hurlbertlab/gdicecco/NLCD_fragmentation/landcover-finescale/")
#codes0111 <- read.csv("nlcd_2001-2011_landcover_change_codes.csv", stringsAsFactors = F)

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
 # area.df <- rasterToPoints(zones.sub, fun = function(x) {x %in% codes0111$ID})
 # write.csv(area.df, paste0("nlcd_30x30_2001_2006_bcr_", bcr, ".csv"), row.names = F)
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
#  area.df <- rasterToPoints(zones.sub, fun = function(x) {x %in% codes0111$ID})
#  write.csv(area.df, paste0("nlcd_30x30_2006_2011_bcr_", bcr, ".csv"), row.names = F)
#}
