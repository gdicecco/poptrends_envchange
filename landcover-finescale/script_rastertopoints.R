## Extract land cover values as tables from rasters of BCRs

#### Libraries ####
library(raster)
library(stringr)

# 2001-2011
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/")
files <- list.files()
nlcd.files <- files[str_detect(files, ".grd")]
grd.files <- nlcd.files[!str_detect(nlcd.files, ".csv")]
dir.grd <- getwd()
nlcd.names <- strsplit(grd.files, "\\.")

for(i in 9:length(grd.files)){
  file <- grd.files[i]
  name <- nlcd.names[[i]][1]
  
  raster <- raster(paste0(dir.grd, "/", file))
  
  df <- as.data.frame(rasterToPoints(raster))
  colnames(df)[3] <- "code"
  write.csv(df, paste0(name, ".csv"), row.names = F)
  
  print(paste0("Completed file no. ", i, " - ", file))
}
