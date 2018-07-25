library(dplyr)
library(SDMTools)
library(raster)
library(stringr)

setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/01-06/")

files <- list.files(pattern = "2006_bcr_")
files.grd <- grep("grd", files)
wd <- "/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/csvs/"

for(i in files.grd) {
  filename <- files.grd[i]
  nameid <- word(filename, sep = "\\.")[1]
  data <- raster(filename)
  fragstats <- ClassStat(data, cellsize = 30, bkgd = 0)
  write.csv(fragstats, paste0(wd, nameid, "_frag.csv"), row.names = F)
}
