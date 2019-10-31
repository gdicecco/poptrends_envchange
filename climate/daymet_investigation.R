# DAYMET investigation

climate_trends <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)
routeDAYMET <- read.csv("bbs_routes_breeding_season_climate.csv", stringsAsFactors = F)

ca_routes <- read_sf("bbs_canada_route_areas.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry) %>%
  mutate(country = "Canada")
us_routes <- read_sf("bbsroutes_5km_buffer.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry) %>%
  mutate(country = "US")
na_routes <- rbind(ca_routes, us_routes)

us_sp <- readOGR("\\\\Bioark.bio.unc.edu/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/bbsroutes_5km_buffer.shp")

test <- filter(na_routes, rteno %in% c(72067, 72066))
plot(test)

df1 <- filter(routeDAYMET, stateroute == 72067)
df2 <- filter(routeDAYMET, stateroute == 72066)

plot(df1$year, df1$tmin)
points(df2$year, df2$tmin, col = "blue")

plot(df1$year, df1$tmax)
points(df2$year, df2$tmax, col = "blue")

test_latlon <- st_transform(test, 4236)

setwd("C:/Users/gdicecco/Desktop/")

download_daymet_ncss(location = c(41, -78, 40, -77),
                     start = 1990,
                     end = 2017,
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
  
  na_routes_transf <- st_transform(test, "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0")
  
  routeclim <- raster::extract(daymet_mean, na_routes_transf, fun = mean, df = T)
  print(routeclim)
  
}

# Pull DAYMET and PRISM for BBS route points and compare in 2017

library(prism)
library(daymetr)
library(sf)
library(raster)
library(ncdf4)
library(tidyverse)

routes <- read.csv("\\\\BioArk//hurlbertlab//databases//BBS//2017//bbs_routes_20170712.csv")

# Daymet
setwd("C:/Users/gdicecco/Desktop/")

# download_daymet_ncss(location = c(60, -145, 15, -52),
#                      start = 2017,
#                      end = 2017,
#                      param = c("tmin", "tmax"), 
#                      frequency = "monthly",
#                      path = "daymet/")

# Read in data
files <- list.files("daymet/")

daymet_nc <- nc_open(paste0("daymet/", files[1]))
daymet_raster <- brick(paste0("daymet/", files[1]))

daymet_may <- daymet_raster[[6]]

daymet_bbs <- raster::extract(daymet_may, dplyr::select(routes, longitude, latitude), df = T)

daymet_routes <- daymet_bbs %>% bind_cols(routes)

names(daymet_routes)[2] <- "daymet_tmax"

# PRISM

options(prism.path = "C:/Users/gdicecco/Desktop/prism/")
get_prism_monthlys("tmax", years = 2017, mon = c(6), keepZip = F)
prism_files <- ls_prism_data()

prism <- prism_stack(prism_files$files[[1]])

prism_bbs <- raster::extract(prism, dplyr::select(daymet_routes, longitude, latitude), df = T) %>%
  bind_cols(daymet_routes)

names(prism_bbs)[2] <- "prism_tmax"

plot(prism_bbs$prism_tmax, prism_bbs$daymet_tmax)
abline(a = 0, b = 1)

cor(prism_bbs$prism_tmax, prism_bbs$daymet_tmax, use = "pairwise.complete.obs")