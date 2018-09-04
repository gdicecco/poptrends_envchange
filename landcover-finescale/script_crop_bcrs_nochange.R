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
setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_1992_landcover_2018_08_31/")
region <- raster("nlcd_1992_30meter_whole.img")

crs.nlcd <- crs(region)
bcrs.proj <- spTransform(us.proj, crs.nlcd)

# output: raster file of landcover change
# output: dataframe of landcover change of interest
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_92/")
for(i in 1:length(bcrs)) {
  bcr <- bcrs[i]
  bcr.sub <- us.proj[us.proj@data$BCR == bcr, ]
  area.sub <- crop(region, extent(bcr.sub))
  zones.sub <- extract(area.sub, bcr.sub, buffer = 200)
  filename <- paste0("nlcd_30x30_1992_2001_bcr_", bcr, ".grd")
  writeRaster(zones.sub, filename = filename)
}
