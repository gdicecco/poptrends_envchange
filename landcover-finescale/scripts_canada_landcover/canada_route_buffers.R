setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/Databases/BBS/GPS_stoplocations/")

library(sf)
library(tidyverse)
library(rgeos)

ca_routes <- read_sf("BBS_routes_USandCanada/BBS_routes_USandCanada.shp") %>%
  filter(COUNTRY == "Canada")

routes <- read.csv("\\\\Bioark.bio.unc.edu/hurlbertlab/Databases/BBS/2017/bbs_routes_20170712.csv") %>%
  mutate(stateroute = statenum*1000 + route) %>%
  dplyr::select(stateroute, latitude, longitude)

ca_lengths <- ca_routes %>%
  mutate(route_length = st_length(.)) %>%
  filter(route_length <= set_units(42000, m) & route_length >= set_units(38000, m)) %>%
  left_join(routes, by = c("rteno" = "stateroute")) %>%
  filter(latitude <= 60)
# 888 routes that meet length and latitude filters

bufferRoutes <- st_buffer(ca_lengths, dist = 5000) %>%
  mutate(route_area = st_area(.))

st_write(bufferRoutes, "C:/Users/gdicecco/Desktop/data/bbs_canada_route_areas.shp")
