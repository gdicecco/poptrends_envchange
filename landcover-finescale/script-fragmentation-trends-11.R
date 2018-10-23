library(dplyr)
library(SDMTools)
library(raster)
library(stringr)

# 2011
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/bcrs_11/")

files <- list.files(pattern = "2011_bcr_")
files.grd <- files[grepl("grd", files)]
wd <- "/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles_nochange/csvs/"

for(i in c(8:9)) {
  filename <- files.grd[i]
  nameid <- word(filename, sep = "\\.")[1]
  data <- raster(filename)
  fragstats <- ClassStat(data, cellsize = 30, bkgd = NA)
  write.csv(fragstats, paste0(wd, nameid, "_frag.csv"), row.names = F)
}
