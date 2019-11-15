## Get climate change data for routes

# 1990 to 2017
# Average monthly min and max temperature, total precipitation during breeding season May - July
# Slope of OLS line = amount of trend

# For routes being used, extract monthly averages, calculate trend

library(sf)
library(rgdal)
library(rgeos)
library(daymetr)
library(tidyverse)
library(raster)
library(ncdf4)
library(lubridate)
library(units)
library(purrr)

## Get BBS routes

## Get route buffers

setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")

ca_routes <- read_sf("bbs_canada_route_areas.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry) %>%
  mutate(country = "Canada")
us_routes <- read_sf("bbsroutes_5km_buffer.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry) %>%
  mutate(country = "US")

na_routes <- rbind(ca_routes, us_routes)

## Get daymet data for north america

# Function: get average temperature for one BBS route polygon

daymetMean <- function(stateroute) {
  rte <- filter(na_routes_transf, rteno == stateroute)
  
  daymet_crop <- raster::crop(daymet_mean, rte)
  daymet_mask <- raster::mask(daymet_crop, rte)
  
  local_extract <- raster::extract(daymet_mask, rte, fun = mean, na.rm = T)
  
  return(local_extract)
}

possibly_daymetMean <- possibly(daymetMean, otherwise = NA)

setwd("C:/Users/gdicecco/Desktop/")

years <- c(1990:2017)

for(y in years) {
  download_daymet_ncss(location = c(60, -145, 15, -52),
                       start = y,
                       end = y,
                       param = c("tmin", "tmax"), 
                       frequency = "monthly",
                       path = "daymet/")
  
  # Read in data
  files <- list.files("daymet/")
  
  for(f in files) {
    daymet_nc <- nc_open(paste0("daymet/", f))
    daymet_raster <- brick(paste0("daymet/", f))
    crs(daymet_raster) <- "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0"
    
    daymet_breeding <- daymet_raster[[5:7]]
    
    daymet_mean <- mean(daymet_breeding, na.rm = T)
    
    crs_daymet <- crs(daymet_breeding)
    
    na_routes_transf <- st_transform(na_routes, "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0")
    
    routeclim <- data.frame(stateroute = na_routes_transf$rteno) %>%
      mutate(mean_temp = purrr::map(stateroute, possibly_daymetMean))
    
    write.csv(routeclim, paste0("C:/Users/gdicecco/Desktop/daymet_out/", f, ".csv"), row.names = F)
    print(f)
    nc_close(daymet_nc)
  }
  
  print(y)
  sapply(paste0("daymet/", files), unlink)

}

# Read in indiv year files
# Join together

dir <- "C:/Users/gdicecco/Desktop/daymet_out/"

routeDAYMET <- data.frame(stateroute = c(), year = c(), mean_tmax = c(), mean_tmin = c())

for(y in years) {
  files <- list.files(dir)
  files_y <- files[grepl(y, files)]
  
  tmax <- read.csv(paste0(dir, files_y[grepl("tmax", files_y)])) %>%
    group_by(stateroute) %>%
    summarize(mean_tmax = mean(mean_temp, na.rm = T))
  tmin <- read.csv(paste0(dir, files_y[grepl("tmin", files_y)])) %>%
    group_by(stateroute) %>%
    summarize(mean_tmin = mean(mean_temp, na.rm = T))
  
  tmp <- tmax %>%
    left_join(tmin) %>%
    mutate(year = y)
  print(nrow(tmp))
  
  routeDAYMET <- rbind(routeDAYMET, tmp)
}

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
# write.csv(routeDAYMET, "bbs_routes_breeding_season_climate.csv", row.names = F)

routeDAYMET <- read.csv("bbs_routes_breeding_season_climate.csv", stringsAsFactors = F)

# Confirm DAYMET means match PRISM
routePRISM <- read.csv("bbs_routes_breeding_season_climate_prism.csv", stringsAsFactors = F)

compare_climate <- routeDAYMET %>%
  left_join(routePRISM, by = c("stateroute", "year"))

cor(compare_climate$mean_tmax, compare_climate$tmax, use = "pairwise.complete.obs")
cor(compare_climate$mean_tmin, compare_climate$tmin, use = "pairwise.complete.obs")

# Calculate climate trend at each route
library(broom)
library(purrr)

route_countries <- na_routes %>%
  dplyr::select(rteno, country) %>%
  st_set_geometry(NULL)

climate_trends <- routeDAYMET %>%
  gather(key = "env", value = "val", c("tmax", "tmin")) %>%
  dplyr::select(-ID) %>%
  left_join(route_countries, by = c("stateroute" = "rteno")) %>%
  group_by(country, stateroute, env, year) %>%
  summarize(breedingAvg = mean(val)) %>%
  group_by(stateroute, env) %>%
  nest() %>%
  mutate(nRow = purrr::map_dbl(data, ~{
    df <- .
    nrow(na.omit(df))
  })) %>%
  filter(nRow != 0) %>%
  mutate(lmFit = purrr::map(data, ~{
    df <- .
    country <- unique(df$country)
    if(country == "Canada") {
      lm(breedingAvg ~ year, filter(df, year >= 1990 & year <= 2010))
    } else {
      lm(breedingAvg ~ year, filter(df, year >= 1992 & year <= 2016))
    }
    
  })) %>%
  mutate(lm_broom = purrr::map(lmFit, tidy)) %>%
  mutate(climateTrend = map_dbl(lm_broom, ~{
    df <- .
    slope <- df %>%
    filter(term == "year") %>%
      dplyr::select(estimate)
    slope[[1]]
  })) %>%
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    p <- df %>%
      filter(term == "year") %>%
      dplyr::select(p.value)
    p[[1]]
  }))

climate_trends_write <- climate_trends %>%
  dplyr::select(stateroute, env, climateTrend, trendPval)

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
#write.csv(climate_trends_write, "bbs_routes_climate_trends.csv", row.names = F)
climate_trends <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)
