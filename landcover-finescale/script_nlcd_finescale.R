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
get.file <- function(x) {
  files2 <- list.files(paste0(dir, "/", area.files[x], ""))
  file.path <- files2[str_detect(files2, "area")]
  return(list(folder = area.files[x], file.name = file.path))
}

# Function to merge two raster files with different origins and extents, same resolution and crs
merge.areas <- function(x, y) {
  extent1 <- extent(x)
  extent2 <- extent(y)
  z <- raster(xmn = min(extent1[1], extent2[1]), xmx = max(extent1[2], extent2[2]), # make blank raster file with extent of the two files added
              ymn = min(extent1[3], extent2[3]), ymx = max(extent1[4], extent2[4]), 
              resolution = res(x), crs = crs(x))
  merge1 <- mosaic(z, x, fun = max, tolerance = max(abs(origin(z) - origin(x)))) # function: where rasters overlap keep larger value
  merge2 <- mosaic(merge1, y, fun = max, tolerance = max(abs(origin(merge1) - origin(y)))) # tolerance allows for different origins 
  return(merge2) # returns the two rasters merged into one
}

#### 1992 - 2001 ####

## Crop raster data for each BCR
region <- raster("nlcd_1992_2001_changepixels_30m_30m.grd")
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/landcover-finescale/")
us.proj <- readOGR("BCRS_contiguous_us.shp")

crs.nlcd <- crs(region)
bcrs.proj <- spTransform(us.proj, crs.nlcd)

bcrs <- c(9, 12, 13, 14, 18, 19, 23, 27, 29)

for(i in 1:length(bcrs)) {
  bcr <- bcrs[i]
  bcr.sub <- bcrs.proj[bcrs.proj@data$BCR == bcr, ]
  area.sub <- crop(region, extent(bcr.sub))
  zones.sub <- mask(area.sub, bcr.sub)
  # where to write these files to
  filename <- paste0("nlcd_30x30_1992_2001_bcr_", bcr, ".grd", sep = "")
  writeRaster(zones.sub, filename = filename)
}

