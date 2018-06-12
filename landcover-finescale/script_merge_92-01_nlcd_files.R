## Merge fine scale zone shapefiles - NLCD 1992-2001

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

## Get codes for landcover change of interest
setwd("/proj/hurlbertlab/nlcd_landcover_change/nlcd_1992_to_2001_landcover_change/")
codes <- read.csv("anderson_land_cover_codes.csv", stringsAsFactors = F) # NLCD land cover class codes

fragcodes <- codes %>%
  filter(grepl("Forest|Grassland", anderson))

## Create map of US
files <- list.files()
area.files <- files[str_detect(files, "area")]
dir <- getwd()

file1 <- get.file(1)
data <- raster(paste0(dir, "/", file1$folder, "/", file1$file.name, sep = ""))
data.frag <- data[data %in% fragcodes$modified]

file2 <- get.file(2)
data2 <- raster(paste0(dir, "/", file2$folder, "/", file2$file.name, sep = ""))
data2.frag <- data2[data2 %in% fragcodes$modified]

region <- merge.areas(data.frag, data2.frag)
crs <- crs(region)

for(i in 3:length(area.files)) { 
  if(i == 8) { # area 8 - Michigan - has different projection
    file3 <- get.file(i)
    data3 <- raster(paste0(dir, "/", file3$folder, "/", file3$file.name, sep = ""))
    data3.frag <- data3[data3 %in% fragcodes$modified]
    data3.proj <- projectRaster(data3.frag, crs = crs)
    region <- merge.areas(region, data3.proj)
  } else {
    file3 <- get.file(i)
    data3 <- raster(paste0(dir, "/", file3$folder, "/", file3$file.name, sep = ""))
    data3.frag <- data3[data3 %in% fragcodes$modified]
    region <- merge.areas(region, data3.frag)
  }
}

setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/")
writeRaster(region, "nlcd_1992_2001_changepixels_30m_30m.grd")
