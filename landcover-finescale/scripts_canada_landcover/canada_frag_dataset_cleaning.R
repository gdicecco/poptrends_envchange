# Canada fragmentation dataset cleaning

library(tidyverse)
library(sf)
library(tmap)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
canada <- read.csv("fragmentation_indices_canada.csv", stringsAsFactors = F) # 870 routes

route_zones <- canada %>% 
  dplyr::select(year, stateroute, file) %>% 
  group_by(year, stateroute) %>% 
  summarize(n_zones = as.factor(n_distinct(file)))

zone_duplicates <- route_zones %>%
  filter(n_zones == 2) # 80

setwd("\\\\BioArk\\hurlbertlab\\databases\\bbs\\2017\\")
setwd("/Volumes/hurlbertlab/databases/bbs/2017/")
routes <- read.csv("bbs_routes_20170712.csv", stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route) %>%
  dplyr::select(stateroute, latitude, longitude)

route_zones_sf <- route_zones %>%
  left_join(routes) %>%
  st_as_sf(coords = c("longitude", "latitude"))

tm_shape(route_zones_sf) + tm_dots(col = "n_zones", size = 0.1, palette = c("blue", "red")) + tm_facets(by = "year")

route_zones_1_sf <- route_zones %>%
  filter(n_zones == 1) %>%
  left_join(routes) %>%
  st_as_sf(coords = c("longitude", "latitude"))

tm_shape(route_zones_1_sf) + tm_dots(col = "n_zones", size = 0.1, palette = "black", legend.show = F) + tm_facets(by = "year")
