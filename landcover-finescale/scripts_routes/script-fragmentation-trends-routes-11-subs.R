# Calculate fragmentation measures for routes, subset of BCRS
# 2011

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
setwd("/proj/hurlbertlab/nlcd_landcover/nlcd_2011_landcover_2011_edition_2014_10_10/")
nlcd <- raster("nlcd_2011_landcover_2011_edition_2014_10_10.img")
routes <- readOGR("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BBS_routepaths/bbsroutes_5km_buffer_bcrsub.shp", 
                  p4s = crs(nlcd))

routenos <- routes@data[ , 1]

setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_routes/")
for(i in 1:nrow(routes@data)) {
  rte <- subset(routes, rteno == routenos[i])
  rtenum <- routenos[i]
  nlcd_crop <- crop(nlcd, rte)
  nlcd_mask <- mask(nlcd_crop, rte)
  class <- ClassStat(nlcd_mask)
  filename <- paste0("classStat_nlcd_30x30_2011_route_", rte, ".csv")
  write.csv(class, filename, row.names = F)
}
