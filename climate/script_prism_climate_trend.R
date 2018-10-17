## Get climate change data for routes

# 1990 to 2017
# Average monthly min and max temperature, total precipitation during breeding season May - Sept
# Slope of OLS line = amount of trend

# For routes being used, extract monthly averages, calculate trend

require(raster)
require(tidyverse)
library(prism)
library(rgdal)
library(rgeos)

options(prism.path = "C:/Users/gdicecco/Desktop/data/prism_1990_2018/")

# Get the PRISM data (1 time only)
# Annual precip, monthly temperature mins and maxes for breeding season
get_prism_annual("ppt", years = 1990:2017, keepZip = F)

get_prism_monthlys("tmin", years = 1990:2018, mon = c(5:9), keepZip = F)

get_prism_monthlys("tmax", years = 1990:2018, mon = c(5:9), keepZip = F)

# Extract breeding season average temps 
prism_files <- ls_prism_data()
prism_files$env <- word(prism_files$files, start = 2, sep = "_")
prism_files$date <- word(prism_files$files, start = 5, sep = "_")
prism_files$year <- substr(prism_files$date, 1, 4)
prism_files$month <- substr(prism_files$date, 5, 6)

prism <- prism_stack(filter(prism_files, year != 2018)$files) # leave out provisional 2018 data for now
prismCRS <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

## Get routes 

setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/Databases/BBS/GPS_stoplocations/")

us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)

routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
routes$stateroute <- routes$statenum*1000 + routes$route

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data$rteno > 3000), ]

routes_subs <- routes %>% # subset stateroutes that were filtered by criteria above
  filter(stateroute %in% us_subs@data$rteno)

xy <- dplyr::select(routes_subs, longitude, latitude) # just need lon-lat for extract

spdf <- SpatialPointsDataFrame(coords = xy, data = routes_subs,
                               proj4string = prismCRS) # routes being used

# Extract breeding season climate data at routes
routePRISM <- raster::extract(prism, xy, df = T) %>%
  bind_cols(xy) %>%
  gather(key = "prismFile", value = "val", 2:309)

routePRISM %<>% mutate(env = word(prismFile, start = 2, sep = "_")) %>%
  mutate(date = word(prismFile, start = 5, sep = "_")) %>%
  mutate(year = substr(date, 1, 4)) %>%
  mutate(month = substr(date, 5, 6))

setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\data\\")
write.csv(routePRISM, "bbs_routes_breeding_season_climate.csv", row.names = F)

# Calculate climate trend at each route
library(broom)
library(purrr)

climate_trends <- routePRISM %>%
  left_join(select(routes_subs, latitude, longitude, bcr, stateroute), by = c("longitude" = "longitude", "latitude" = "latitude")) %>%
  group_by(stateroute, env, year) %>%
  mutate(breedingAvg = mean(val)) %>%
  group_by(stateroute, env) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    df.short <- df %>%
      select(year, breedingAvg) %>%
      unique()
    lm(year ~ breedingAvg, df.short)
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  mutate(climateTrend = map_dbl(lm_broom, ~{
    df <- .
    slope <- df %>%
    filter(term == "breedingAvg") %>%
      select(estimate)
    slope[[1]]
  })) %>%
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    p <- df %>%
      filter(term == "breedingAvg") %>%
      select(p.value)
    p[[1]]
  }))

climate_trends_write <- climate_trends %>%
  select(stateroute, env, climateTrend, trendPval)

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
write.csv(climate_trends_write, "bbs_routes_climate_trends.csv", row.names = F)
