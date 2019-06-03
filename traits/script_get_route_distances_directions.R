## Get distance of BBS routes to range center for forest breeding bird species

library(tidyverse)
library(sf)
library(rgdal)
library(geosphere)

## List of species observed in BCRs of interest during time window (1990-present)
routes <- read.csv("\\\\BioArk\\HurlbertLab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

routes <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_routes_20170712.csv")
weather <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_weather_20170712.csv")
counts <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_counts_20170712.csv")
species <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_species_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# species four letter codes
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
spp_codes <- read.csv("four_letter_codes_birdspp.csv", stringsAsFactors = F)
spp_codes$common_name_lower <- tolower(spp_codes$COMMONNAME)

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)
setwd("\\\\BioArk/hurlbertlab/Databases/BBS/GPS_stoplocations/")
setwd("/Volumes/hurlbertlab/Databases/BBS/GPS_stoplocations/")
us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data$rteno > 3000), ]

# Subset stateroutes: runtype = 1, sampled once every 5 years since 1990, between 38 and 42 km in length
routes.short <- RT1.routes %>%
  filter(stateroute %in% us_subs@data$rteno) %>%
  mutate(year_bin = 5*floor(year/5) + 5/2) %>%
  filter(year >= 1990) %>%
  group_by(stateroute) %>%
  distinct(year_bin) %>%
  count() %>%
  filter(n >= 5)
# 1513 routes

# Subset species: forest breeding birds by area sensitivity

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
area_sensitive <- read.csv("traits/area-sensitivity-forest-birds-Bouliner1998.csv", stringsAsFactors = F)
area_sensitive <- area_sensitive[, -c(4)]

species$common_name_lower <- tolower(species$english_common_name)
area_sensitive$common_name_lower <- tolower(area_sensitive$Common_name)

# Species that don't match: Northern Flicker (4120), Chuck-will's-widow (4160), Whip-poor-will (4171), Black-and-white warbler (6360)
spp_aous <- data.frame(common_name_lower = c("northern flicker", "chuck-will’s-widow", "whip-poor-will"),
                       aou = c(4120, 4160, 4171)) %>%
  left_join(area_sensitive) %>%
  left_join(dplyr::select(species, -common_name_lower), by = "aou")

area_aous <- area_sensitive %>%
  left_join(species) %>%
  filter(!is.na(aou)) %>%
  bind_rows(spp_aous) %>%
  left_join(spp_codes) %>%
  filter(!is.na(Wintering_Slimit_general))

fourletter_codes <- dplyr::select(area_aous, aou, SPEC)

## BBS subset

counts.subs <- counts %>%
  filter(aou %in% area_aous$aou) %>%
  merge(filter(RT1.routes, stateroute %in% routes.short$stateroute), by = c("stateroute", "year")) %>%
  filter(year >= 1990, year < 2017)
# 1513 routes

### Join species list with range maps name

rangemap_list <- data.frame(files = list.files(path = "\\\\BioArk/HurlbertLab/GIS/birds/All/All/")) %>%
  mutate(genus = word(files, start = 1, sep = "_"),
         species = word(files, start = 2, sep = "_"),
         id = word(files, start = 3, sep = "_"),
         BL_sciname = paste(genus, species)) %>%
  filter(grepl("shp", files)) %>%
  dplyr::select(files, BL_sciname)

# Try different scientific names to match up species lists
range_spp_list <- area_aous %>%
  left_join(rangemap_list, by = c("Scientific_name" = "BL_sciname")) 

spp_matches <- range_spp_list %>%
  filter(!(is.na(files))) %>%
  dplyr::select(Common_name, aou, SPEC, files)

spp_nas <- range_spp_list %>%
  filter(is.na(files)) %>%
  dplyr::select(-files) %>%
  left_join(rangemap_list, by = c("SCINAME" = "BL_sciname"))

spp_matches_2 <- spp_nas %>%
  filter(!(is.na(files))) %>%
  dplyr::select(Common_name, aou, SPEC, files)

# Species that still don't match up with Bird Life names - fix manually
BL_names <- c("Parus atricapillus", 
              "Parus carolinensis", 
              "Leuconotopicus villosus",
              "Hylatomus pileatus")

spp_still_nas <- spp_nas %>%
  filter(is.na(files)) %>%
  dplyr::select(-files) %>%
  mutate(BL_names = BL_names) %>%
  left_join(rangemap_list, by = c("BL_names" = "BL_sciname")) %>%
  dplyr::select(Common_name, aou, SPEC, files)

spp_map_files <- bind_rows(spp_matches, spp_matches_2, spp_still_nas)

### Get range centroid for each species

path <- c("\\\\BioArk/HurlbertLab/GIS/birds/All/All/")
range_centroids <- data.frame(aou = c(), filename = c(), 
                              centroid_lon = c(), centroid_lat = c())

for(spec in spp_map_files$aou) {
  file <- spp_map_files[spp_map_files$aou == spec, 4]
  range <- read_sf(paste0(path, file)) %>%
    filter(SEASONAL == 1 | SEASONAL == 2 | SEASONAL == 5) %>%
    st_transform('+proj=longlat +ellps=GRS80 +no_defs') %>%
    st_union()
  centroid <- st_coordinates(st_centroid(range))
  
  temp <- data.frame(aou = spec,
                     filename = file,
                     centroid_lon = centroid[1],
                     centroid_lat = centroid[2])
  
  range_centroids <- bind_rows(range_centroids, temp)
}
write.csv(range_centroids, "traits/species_Brange_centroids.csv", row.names = F)

### Distance and direction of each route from range centroid

route_directions <- counts.subs %>%
  filter(aou %in% spp_map_files$aou) %>%
  group_by(aou) %>%
  left_join(range_centroids) %>%
  nest() %>%
  mutate(routes = map(data, ~{
    df <- .
    stateroutes <- unique(df$stateroute)
    
    c_lon <- unique(df$centroid_lon)
    c_lat <- unique(df$centroid_lat)
    
    cent <- c(c_lon, c_lat)
    
    coords <- routes %>%
      filter(stateroute %in% stateroutes) %>%
      dplyr::select(stateroute, latitude, longitude)
    
    results <- data.frame(stateroute = c(), dist_from_centroid = c(), bearing = c())
    
    for(route in coords$stateroute) {
      rt <- filter(coords, stateroute == route)
      pt_route <- c(rt$longitude, rt$latitude)
      
      dist <- distGeo(cent, pt_route)
      bearing <- bearing(cent, pt_route)
      
      res <- data.frame(stateroute = route, dist_from_centroid = dist, bearing = bearing)
      
      results <- bind_rows(results, res)
      
    }
    
    results
    
  }))

route_directions_unnest <- route_directions %>%
  dplyr::select(-data) %>%
  unnest()

write.csv(route_directions_unnest, "traits/species_route_dist_from_centroids.csv",
          row.names = F)
