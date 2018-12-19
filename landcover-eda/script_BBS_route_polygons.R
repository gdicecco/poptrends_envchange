library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(dplyr)
library(stringr)

setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/Databases/BBS/GPS_stoplocations/")

us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)

routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
routes$stateroute <- routes$statenum*1000 + routes$route

routes_subs <- filter(routes, bcr %in% bcrs)

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data$rteno > 3000), ]

us_short <- subset(us_subs, rteno %in% routes_subs$stateroute)

bufferRoutes <- gBuffer(us_subs, width = 5000, byid = TRUE)

setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\")
writeOGR(bufferRoutes, ".", "bbsroutes_5km_buffer", driver = "ESRI Shapefile")

setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")
setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\nlcd_frag_proj_shapefiles\\")
bufferRoutes <- readOGR("bbsroutes_5km_buffer.shp")

library(SDMTools)

setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\GIS\\LandCoverData\\")
setwd("/Volumes/hurlbertlab/GIS/LandCoverData/")
nlcd2001 <- raster("nlcd_2001_landcover_2011_edition_2014_10_10/nlcd_2001_landcover_2011_edition_2014_10_10.img")

bufferRoutes_transf <- spTransform(bufferRoutes, crs(nlcd2001))
bufferRoutes_tBCR <- subset(bufferRoutes_transf, rteno %in% routes_subs$stateroute)

setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\")
writeOGR(bufferRoutes_tBCR, ".", "bbsroutes_5km_buffer_bcrsub", driver = "ESRI Shapefile")

route1 <- subset(bufferRoutes_transf, rteno == 2001)

nlcd2001_subs <- crop(nlcd2001, extent(route1))
route_nlcd <- mask(nlcd2001_subs, route1)

plot(route_nlcd)
