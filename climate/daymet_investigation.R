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

crs(daymet_raster) <- "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0"


daymet_june <- daymet_raster[[6]]
daymet_january <- daymet_raster[[1]]

routes_sf <- dplyr::select(routes, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4236) %>%
  st_transform("+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0")

daymet_bbs <- raster::extract(daymet_june, st_coordinates(routes_sf), df = T)
daymet_bbs_jan <- raster::extract(daymet_january, st_coordinates(routes_sf), df = T)

daymet_routes <- daymet_bbs %>% 
  bind_cols(routes) %>%
  left_join(daymet_bbs_jan)

names(daymet_routes)[c(2,15)] <- c("daymet_tmax_jun", "daymet_tmax_jan")

# PRISM

options(prism.path = "C:/Users/gdicecco/Desktop/prism/")
get_prism_monthlys("tmax", years = 2017, mon = c(1,6), keepZip = F)
prism_files <- ls_prism_data()

routes_points <- dplyr::select(routes, route, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"))

prism <- prism_stack(prism_files$files[c(1,3)])

prism_bbs <- raster::extract(prism, dplyr::select(daymet_routes, longitude, latitude), df = T) 

names(prism_bbs)[2:3] <- c("prism_tmax_jan", "prism_tmax_jun")

tmax_compare <- daymet_routes %>%
  left_join(prism_bbs)

jan_cor <- cor(tmax_compare$prism_tmax_jan, tmax_compare$daymet_tmax_jan, use = "pairwise.complete.obs")
jun_cor <- cor(tmax_compare$prism_tmax_jun, tmax_compare$daymet_tmax_jun, use = "pairwise.complete.obs")

jan_plot <- ggplot(tmax_compare, aes(x = tmax_compare$daymet_tmax_jan, y = tmax_compare$prism_tmax_jan)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) + theme_classic(base_size = 12)
