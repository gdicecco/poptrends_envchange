library(dplyr)
library(SDMTools)
library(raster)
library(stringr)

# 1992
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_92_cloud/")

files <- list.files(pattern = "1992_bcr_")
files.grd <- files[grepl("grd", files)]
wd <- "/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_92_cloud/"

for(i in length(files.grd) {
  filename <- files.grd[i]
  nameid <- word(filename, sep = "\\.")[1]
  data <- raster(filename)
  fragstats <- ClassStat(data, cellsize = 30, bkgd = NA)
  write.csv(fragstats, paste0(wd, nameid, "_frag.csv"), row.names = F)
}
