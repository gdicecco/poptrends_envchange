## Make figures for cartoon methods

# Rasters of tmin, tmax, ppt for one route

require(raster)
require(tidyverse)
library(prism)
library(rgdal)
library(rgeos)

options(prism.path = "/Users/gracedicecco/Desktop/prism_2018/")

# Get the PRISM data (1 time only)
# Annual precip, monthly temperature mins and maxes for breeding season
get_prism_annual("ppt", years = 2017, keepZip = F)

get_prism_monthlys("tmin", years = 2017, mon = c(5), keepZip = F)

get_prism_monthlys("tmax", years = 2017, mon = c(5), keepZip = F)

# Extract breeding season average temps 
prism_files <- ls_prism_data()
prism_files$env <- word(prism_files$files, start = 2, sep = "_")
prism_files$date <- word(prism_files$files, start = 5, sep = "_")
prism_files$year <- substr(prism_files$date, 1, 4)
prism_files$month <- substr(prism_files$date, 5, 6)

prism <- prism_stack(filter(prism_files, year != 2018)$files) # leave out provisional 2018 data for now
prismCRS <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")
bufferRoutes <- readOGR("bbsroutes_5km_buffer.shp")
bufferRoutes_transf <- spTransform(bufferRoutes, prismCRS)

route1 <- subset(bufferRoutes_transf, rteno == 2001)

list2env(setNames(unstack(prism), names(prism)), .GlobalEnv)


library(tmap)
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/figures/")
prism_crop <- crop(PRISM_ppt_stable_4kmM3_2017_bil, route1)
route_prism <- mask(prism_crop, route1)
tm_shape(route_prism) +
  tm_raster(style = "cont", title = "Ppt", palette = "YlGn")
ggsave("methods_ppt_raster.tiff", units = "in")

prism_crop_tmax <- crop(PRISM_tmax_stable_4kmM2_201705_bil, route1)
route_prism_tmax <- mask(prism_crop_tmax, route1)
tm_shape(route_prism_tmax) +
  tm_raster(style = "cont", title = "Tmax", palette = "Reds", legend.show = T)

prism_crop_tmin <- crop(PRISM_tmin_stable_4kmM2_201705_bil, route1)
route_prism_tmin <- mask(prism_crop_tmin, route1)
tm_shape(route_prism_tmax) +
  tm_raster(style = "cont", title = "Tmin", palette = "Blues")

# Plots of abundance trend and climate trend for one species

setwd("/Volumes/hurlbertlab/DiCecco/data/")
routePRISM <- read.csv("bbs_routes_breeding_season_climate.csv", stringsAsFactors = F)
routeClim <- filter(routePRISM, ID == 6) %>%
  group_by(year, env) %>%
  summarize(mean = mean(val))

library(ggplot2)
theme_set(theme_classic())

ggplot(filter(routeClim, env == "ppt"), aes(x = year, y = mean)) + geom_point() + 
  geom_smooth(method = "lm", se = F, color = "green") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(x = "Year", y = "Precipitation (mm)")

ggplot(filter(routeClim, env == "tmax"), aes(x = year, y = mean)) + geom_point() + 
  geom_smooth(method = "lm", se = F, color = "red") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(x = "Year", y = "Tmax (C)")

ggplot(filter(routeClim, env == "tmin"), aes(x = year, y = mean)) + geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(x = "Year", y = "Tmin (C)")

routes <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_routes_20170712.csv")
counts <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_counts_20170712.csv")
species <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_species_20170712.csv")
weather <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_weather_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

counts_onespp <- counts %>%
  filter(aou == 2890, year > 1989, stateroute == 2001)

ggplot(counts_onespp, aes(x = year, y = speciestotal)) + geom_point() + 
  geom_smooth(method = "lm", se = F, color = "purple") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(x = "Year", y = "Abundance")
