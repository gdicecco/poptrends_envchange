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
library(sf)

options(prism.path = "C:/Users/gdicecco/Desktop/data/prism_1990_2018/")

# Get the PRISM data (1 time only)
# Annual precip, monthly temperature mins and maxes for breeding season

get_prism_monthlys("tmin", years = 1990:2017, mon = c(5:7), keepZip = F)

get_prism_monthlys("tmax", years = 1990:2017, mon = c(5:7), keepZip = F)

# Extract breeding season average temps 
prism_files <- ls_prism_data()
prism_files$env <- word(prism_files$files, start = 2, sep = "_")
prism_files$date <- word(prism_files$files, start = 5, sep = "_")
prism_files$year <- substr(prism_files$date, 1, 4)
prism_files$month <- substr(prism_files$date, 5, 6)

## Get routes 

setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")
  
us_routes <- read_sf("bbsroutes_5km_buffer.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry) %>%
  mutate(country = "US")

na_routes <- us_routes %>%
  st_transform("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Extract breeding season climate data at routes (group by year and average)
years <- c(1990:2017)

prism_df <- data.frame(stateroute = c(),
                       year = c(),
                       tmax = c(),
                       tmin = c())

for(y in years) {
  prism <- prism_stack(filter(prism_files, year == y, grepl("tmax", env))$files)
  prism_mean <- mean(prism, na.rm = T)
  routePRISM <- raster::extract(prism_mean, na_routes, fun = mean, na.rm = T, df = T) 
  
  prism_tmin <- prism_stack(filter(prism_files, year == y, grepl("tmin", env))$files)
  prism_tmin_mean <- mean(prism_tmin, na.rm = T)
  routePRISM_tmin <- raster::extract(prism_tmin_mean, na_routes, fun = mean, na.rm = T, df = T)
  
  prism_df <- rbind(prism_df, data.frame(stateroute = na_routes$rteno,
                                         year = y,
                                         tmax = routePRISM$layer,
                                         tmin = routePRISM_tmin$layer))
  print(y)
  
}

write.csv(prism_df, "climate/bbs_routes_breeding_season_climate_prism.csv", row.names = F)

routePRISM <- read.csv("climate/bbs_routes_breeding_season_climate_prism.csv", stringsAsFactors = F)

# Calculate climate trend at each route
library(broom)
library(purrr)

climate_trends <- routePRISM %>%
  gather(key = env, value = breedingAvg, 3:4) %>%
  group_by(stateroute, env) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    df.short <- df %>%
      filter(year >= 1992 & year <= 2016) %>%
      select(year, breedingAvg) %>%
      unique()
    lm(breedingAvg ~ year, df.short)
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  mutate(climateTrend = map_dbl(lm_broom, ~{
    df <- .
    slope <- df %>%
      filter(term == "year") %>%
      select(estimate)
    slope[[1]]
  })) %>%
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    p <- df %>%
      filter(term == "year") %>%
      select(p.value)
    p[[1]]
  }))

climate_trends_write <- climate_trends %>%
  select(stateroute, env, climateTrend, trendPval)

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
write.csv(climate_trends_write, "bbs_routes_climate_trends_prism.csv", row.names = F)
climate_trends <- read.csv("bbs_routes_climate_trends_prism.csv", stringsAsFactors = F)

climate_trends_daymet <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)

climate_trends_compare <- climate_trends %>%
  left_join(climate_trends_daymet, by = c("stateroute", "env"), suffix = c("_prism", "_daymet"))

cor(climate_trends_compare$climateTrend_daymet, climate_trends_compare$climateTrend_prism)
