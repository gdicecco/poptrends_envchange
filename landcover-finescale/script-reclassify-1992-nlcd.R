setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\GIS\\LandCoverData\\nlcd_1992_landcover_2018_08_31\\")
nlcd1992 <- raster("nlcd_1992_whole.tif")

nlcd1992_conv <- reclassify(nlcd1992, c(32,33,31, 61,61,82, 83,84,82))

setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\")
writeRaster(nlcd1992_conv, "nlcd_1992_whole_reclassified.tif")
