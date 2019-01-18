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

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/")
ggplot(route_trends, aes(x = tmax, y = dEDz)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmax", y = "Trend in deltaED")
ggsave("route_tmax_ded.tiff", units = "in")

ggplot(route_trends, aes(x = tmin, y = dEDz)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmin", y = "Trend in deltaED")
ggsave("route_tmin_ded.tiff", units = "in")

ggplot(route_trends, aes(x = tmin, y = ppt)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmin", y = "Trend in ppt")
ggsave("route_tmin_ppt.tiff", units = "in")

ggplot(route_trends, aes(x = tmax, y = ppt)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmax", y = "Trend in ppt")
ggsave("route_tmax_ppt.tiff", units = "in")

ggplot(route_trends, aes(x = tmax, y = tmin)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmax", y = "Trend in Tmin")
ggsave("route_tmax_tmin.tiff", units = "in")


# Map of route changes

library(sf)
setwd("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us_sf <- read_sf("BCRs_contiguous_us.shp")

routes_sf <- st_as_sf(route_trends, coords = c("longitude", "latitude"))

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/")
us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "gray")
tmax_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "tmax", palette = "-RdBu", midpoint = NA, size = 0.2, style = "cont", title = "Tmax")
tmax_map
tmap_save(tmax_map, "routes_tmax_map.tiff", units = "in")

tmin_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "tmin", palette = "-RdBu", midpoint = NA, size = 0.2, style = "cont", title = "Tmin")
tmin_map
tmap_save(tmin_map, "routes_tmin_map.tiff", units = "in")
# flip palette direction

ppt_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "ppt", palette = "PRGn", midpoint = NA, size = 0.2, style = "cont", title = "Ppt")
ppt_map
tmap_save(ppt_map, "routes_ppt_map.tiff", units = "in")

ded_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "dEDz", palette = "-PiYG", midpoint = NA, size = 0.2, style = "cont", title = "deltaED")
ded_map
tmap_save(ded_map, "routes_dED_map.tiff", units = "in")

## Make plots like this for forest fragmentation
head(frags)

classlegend00s <- data.frame(class = c(11:12, 21:24, 31, 41:43, 51, 52, 71, 81:82, 90, 95), 
                             legend = c("Open water", "Perennial ice/snow", "Developed, open space", "Developed, low intensity", "Developed, medium intensity", "Developed, high intensity", "Barren land", "Deciduous forest", "Evergreen forest", "Mixed forest", "Dwarf scrub", "Shrub/scrub", "Grassland/herbaceous", "Pasture/hay", "Cultivated crops", "Woody wetlands", "Emergent herbaceous wetlands"))
classlegend92 <- data.frame(class = c(11:12, 85, 21:23, 31:33, 41:43, 51, 61, 71, 81:84, 91:92),
                            legend = c("Open water", "Perennial ice/snow", "Developed, open space", "Developed, low intensity", "Developed, medium intensity", "Developed, high intensity", "Barren land", "Barren land", "Barren land", "Deciduous forest", "Evergreen forest", "Mixed forest", "Shrub/scrub", "Cultivated crops", "Grassland/herbaceous", "Pasture/hay",  "Cultivated crops",  "Cultivated crops",  "Cultivated crops", "Woody wetlands", "Emergent herbaceous wetlands"))

frags.92 <- frags %>%
  filter(year == 1992) %>%
  left_join(classlegend92)
frags.00s <- frags %>%
  filter(year > 1992) %>%
  left_join(classlegend00s)

forest_ed <- bind_rows(frags.92, frags.00s) %>%
  filter(grepl("forest", legend)) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge)/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  mutate(deltaED = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))

forest_ed$dEDz <- (forest_ed$deltaED - mean(na.omit(forest_ed$deltaED)))/sd(na.omit(forest_ed$deltaED))

forest_map <- us + tm_shape(forest_ed) + 
  tm_dots(col = "dEDz", palette = "-RdYlGn", midpoint = NA, size = 0.2, style = "cont", title = "deltaED Forest")
forest_map 
tmap_save(forest_map, "routes_dED_forest_map.tiff", units = "in")

# Urbanization
urban <- bind_rows(frags.92, frags.00s) %>%
  filter(grepl("Developed", legend)) %>%
  group_by(stateroute, year) %>%
  summarize(urban.prop = sum(prop.landscape)) %>%
  spread(key = "year", value = "urban.prop") %>%
  mutate(deltaUP = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))
urban$dUPz <- (urban$deltaUP - mean(na.omit(urban$deltaUP)))/sd(na.omit(urban$deltaUP))

urban_map <- us + tm_shape(urban) + 
  tm_dots(col = "deltaUP", palette = "YlOrRd", midpoint = NA, size = 0.2, style = "cont", title = "dUrban (proportion of landscape)")
urban_map
tmap_save(urban_map, "routes_urban_map.tiff", units = "in")

# PCA of env. change on routes
head(route_trends)

pca <- prcomp(route_trends[,c(2:4, 6)], scale = T)
summary(pca)
biplot(pca)

              