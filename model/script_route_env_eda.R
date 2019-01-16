## Analysis of route-level environmental change (climate, habitat fragmentation)

library(tidyverse)
library(raster)
library(rgdal)
library(tmap)
library(ggplot2)
library(cowplot)

######## Reading in and subsetting data ##########
# Population data
## BBS
## List of species observed in BCRs of interest during time window (1990-present)
routes <- read.csv("\\\\BioArk\\HurlbertLab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

routes <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_routes_20170712.csv")
weather <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_weather_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)
setwd("\\\\BioArk/hurlbertlab/Databases/BBS/GPS_stoplocations/")
setwd("/Volumes/hurlbertlab/Databases/BBS/GPS_stoplocations/")
us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data$rteno > 3000), ]

# plot of routes

setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
setwd("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us.proj <- readOGR("BCRs_contiguous_us.shp")

us_subs_transf <- spTransform(us_subs, crs(us.proj))

plot(us.proj, col = "gray73", border = "gray73")
plot(us_subs_transf, add = T)

routes.short <- RT1.routes %>% # subset stateroutes that were filtered by criteria above
  filter(stateroute %in% us_subs@data$rteno)
# 2161 routes

# Habitat fragmentation data
frags <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_nlcd.csv", stringsAsFactors = F)

route_ed <- frags %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge)/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2011` - `1992`) %>%
  mutate(dEDz = (deltaED - mean(deltaED))/sd(deltaED))

# Climate data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
climate_trends <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)

climate_wide <- climate_trends %>%
  dplyr::select(-trendPval) %>%
  spread(key = "env", value = "climateTrend")

route_trends <- climate_wide %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(route_ed) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude))

####### Route level trends in environmental change

ggplot(route_trends, aes(x = tmax, y = dEDz)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmax", y = "Trend in deltaED")

# Map of route changes

library(sf)
setwd("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us_sf <- read_sf("BCRs_contiguous_us.shp")

routes_sf <- st_as_sf(route_trends, coords = c("longitude", "latitude"))

us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "gray")
us + tm_shape(routes_sf) + 
  tm_dots(col = "tmax", palette = "RdBu", midpoint = NA, size = 0.2, style = "cont", title = "Tmax")

us + tm_shape(routes_sf) + 
  tm_dots(col = "tmin", palette = "RdBu", midpoint = NA, size = 0.2, style = "cont", title = "Tmin")

us + tm_shape(routes_sf) + 
  tm_dots(col = "ppt", palette = "PRGn", midpoint = NA, size = 0.2, style = "cont", title = "Ppt")

us + tm_shape(routes_sf) + 
  tm_dots(col = "dEDz", palette = "PiYG", midpoint = NA, size = 0.2, style = "cont", title = "deltaED")

## Make plots like this for forest fragmentation and urbanization (?)

# PCA of env. change on routes
head(route_trends)

pca <- prcomp(route_trends[,c(2:4, 6)], scale = T)
summary(pca)
biplot(pca)

              