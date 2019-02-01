library(raster)

setwd("\\\\BioArk\\hurlbertlab\\GIS\\LandCoverData\\")

nlcd1992 <- raster("nlcd_1992_landcover_2018_08_31\\nlcd_1992_whole.tif")

nlcd2001 <- raster("nlcd_2001_landcover_2011_edition_2014_10_10\\nlcd_2001_landcover_2011_edition_2014_10_10.img")

nlcd2006 <- raster("nlcd_2006_landcover_2011_edition_2014_10_10\\nlcd_2006_landcover_2011_edition_2014_10_10.img")

nlcd2011 <- raster("nlcd_2011_landcover_2011_edition_2014_10_10\\nlcd_2011_landcover_2011_edition_2014_10_10.img")

## Reclassify to Urban, Forest, Shrubland, Grassland, Agriculture, Wetlands
newcode <- data.frame(code = seq(1,9), 
                      legend = c("Open water", "Urban", "Barren", "Forest", "Shrubland", "Agricultural", "Grasslands", "Wetlands", "Perennial ice, snow"))

# Reclassify 1992 data
## Create reclassification matrix
reclass_simp_92 <- matrix(nrow = 21, ncol = 2)
reclass_simp_92[, 1] <- c(11,12,85,21:23,31:33,41:43,51,61,71,81:84,91:92) # NLCD 1992 class codes
reclass_simp_92[, 2] <- c(1, 9, rep(2,4), rep(3, 3), rep(4,3), 5, 6, 7, rep(6, 4), rep(8, 2))
colnames(reclass_simp_92) <- c("is", "becomes")

nlcd1992_simp <- reclassify(nlcd1992, rcl = reclass_simp_92, right = NA)

writeRaster(nlcd1992_simp, "nlcd_1992_landcover_2018_08_31\\nlcd_1992_whole_simplified.tif", overwrite = T)

# Reclassify 2001, 2006, 2011 data
## Create reclassification matrix
reclass_simp_00s <- matrix(nrow = 17, ncol = 2)
reclass_simp_00s[, 1] <- c(11, 12, 21:24, 31, 41:43, 52, 82, 71, 81:82, 90, 95)
reclass_simp_00s[, 2] <- c(1, 9, rep(2,4), 3, rep(4, 3), 5, 6, 7, rep(6, 2), rep(8, 2))
colnames(reclass_simp_00s) <- c("is", "becomes")

nlcd2001_simp <- reclassify(nlcd2001, rcl = reclass_simp_00s, right = NA)

writeRaster(nlcd2001_simp, "nlcd_2001_landcover_2011_edition_2014_10_10\\nlcd_2001_whole_simplified.tif", overwrite = T)

nlcd2006_simp <- reclassify(nlcd2006, rcl = reclass_simp_00s, right = NA)

writeRaster(nlcd2006_simp, "nlcd_2006_landcover_2011_edition_2014_10_10\\nlcd_2006_whole_simplified.tif", overwrite = T)

nlcd2011_simp <- reclassify(nlcd2011, rcl = reclass_simp_00s, right = NA)

writeRaster(nlcd2001_simp, "nlcd_2011_landcover_2011_edition_2014_10_10\\nlcd_2011_whole_simplified.tif", overwrite = T)
