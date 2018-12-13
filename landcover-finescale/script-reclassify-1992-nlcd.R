library(raster)

setwd("\\\\BioArk\\hurlbertlab\\GIS\\LandCoverData\\nlcd_1992_landcover_2018_08_31\\")
nlcd1992 <- raster("nlcd_1992_whole.tif")

reclass <- matrix(nrow = 5, ncol = 2)
reclass[, 1] <- c(32, 33, 61, 83, 84)
reclass[, 2] <- c(31, 31, 82, 82, 82)
colnames(reclass) <- c("is", "becomes")

nlcd1992_reclass <- reclassify(nlcd1992, rcl = reclass, right = NA)

writeRaster(nlcd1992_reclass, "nlcd_1992_whole_reclassified.tif", overwrite = T)
