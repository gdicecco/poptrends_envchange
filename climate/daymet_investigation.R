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
