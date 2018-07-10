## Extract land cover values as tables from rasters of BCRs

#### Libraries ####
library(raster)
library(stringr)

# 1992-2001
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/92-01/")
files <- list.files()
nlcd.files <- files[str_detect(files, ".grd")]
dir.grd <- getwd()
nlcd.names <- strsplit(nlcd.files, "\\.")

setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/csvs/")
for(i in 1:length(nlcd.files)) {
  file <- nlcd.files[i]
  name <- nlcd.names[[i]][1]
  
  raster <- raster(paste0(dir.grd, "/", file))
  
  df <- as.data.frame(rasterToPoints(raster))
  colnames(df)[3] <- "code"
  write.csv(df, paste0(name, ".csv"), row.names = F)
  
  print(paste0("Completed file no. ", i, " - ", file))
}

# 2001-2006
setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/01-06/")
files <- list.files()
nlcd.files <- files[str_detect(files, ".grd")]
dir.grd <- getwd()
nlcd.names <- strsplit(nlcd.files, "\\.")

setwd("/proj/hurlbertlab/gdicecco/nlcd_frag_proj_shapefiles/BCRs_landcover_output/csvs/")
for(i in 7:8){
  file <- nlcd.files[i]
  name <- nlcd.names[[i]][1]
  
  raster <- raster(paste0(dir.grd, "/", file))
  
  df <- as.data.frame(rasterToPoints(raster))
  colnames(df)[3] <- "code"
  write.csv(df, paste0(name, ".csv"), row.names = F)
  
  print(paste0("Completed file no. ", i, " - ", file))
}
