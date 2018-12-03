# Calculate fragmentation measures for routes, all BCRS
# 1992

#### Libraries ####
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(dplyr)
library(stringr)
library(SDMTools)

# Read in BBS routes shapefile

# 2001
setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_1992_landcover_2018_08_31/")
nlcd <- raster("nlcd_1992_whole_reclassified.tif")
routes <- readOGR("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BBS_routepaths/bbsroutes_5km_buffer.shp")
routes_tr <- spTransform(routes, crs(nlcd))

routenos <- routes_tr@data[ , 1]

setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_routes_all/")
for(i in 1:nrow(routes_tr@data)) {
  rte <- subset(routes_tr, rteno == routenos[i])
  rtenum <- routenos[i]
  nlcd_crop <- crop(nlcd, rte)
  nlcd_mask <- mask(nlcd_crop, rte)
  class <- ClassStat(nlcd_mask)
  filename <- paste0("classStat_nlcd_30x30_1992_route_", rtenum, ".csv")
  write.csv(class, filename, row.names = F)
}