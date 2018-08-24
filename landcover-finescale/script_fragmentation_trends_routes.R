## Fragmentation trends, routes w/ 500 m buffer only

library(dplyr)
library(SDMTools)
library(raster)
library(stringr)

# 2001-2006
#setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/01-06/")

#files <- list.files(pattern = "2006_routes_")
#files.grd <- files[grepl("grd", files)]
#wd <- "/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/csvs/"

#for(i in 5:6) {
#  filename <- files.grd[i]
#  nameid <- word(filename, sep = "\\.")[1]
#  data <- raster(filename)
#  fragstats <- ClassStat(data, cellsize = 30, bkgd = 0)
#  write.csv(fragstats, paste0(wd, nameid, "_frag.csv"), row.names = F)
#}

# 2006-2011
#setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/06-11/")

#files <- list.files(pattern = "2011_routes_")
#files.grd <- files[grepl("grd", files)]
#wd <- "/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/csvs/"

#for(i in 1:length(files.grd)) {
#  filename <- files.grd[i]
#  nameid <- word(filename, sep = "\\.")[1]
#  data <- raster(filename)
#  fragstats <- ClassStat(data, cellsize = 30, bkgd = 0)
#  write.csv(fragstats, paste0(wd, nameid, "_frag.csv"), row.names = F)
#}

# 1992-2001
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/92-01/")

files <- list.files(pattern = "2001_routes_")
files.grd <- files[grepl("grd", files)]
wd <- "/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/csvs/"

for(i in 1:9) {
  filename <- files.grd[i]
  nameid <- word(filename, sep = "\\.")[1]
  data <- raster(filename)
  fragstats <- ClassStat(data, cellsize = 30, bkgd = NA)
  write.csv(fragstats, paste0(wd, nameid, "_frag.csv"), row.names = F)
}

