
# source script_climate_trends

daymet_1990 <- filter(routeDAYMET, year == 2000) %>% dplyr::select(year, stateroute, mean_tmax)

na_daymet <- na_routes_transf %>%
  left_join(daymet_1990, by = c("rteno" = "stateroute"))

results <- data.frame(stateroute = c(), global_extract = c(), local_extract = c(), local_mean = c())

for(rte in na_daymet$rteno) {
  test_rte <- filter(na_daymet, rteno == rte)
  daymet_crop <- raster::crop(daymet_mean, test_rte)
  daymet_mask <- raster::mask(daymet_crop, test_rte)
  
  local_extract <- raster::extract(daymet_mask, test_rte, fun = mean, na.rm = T)
  
  daymet_df <- rasterToPoints(daymet_mask)
  local_mean <- mean(daymet_df[, 3], na.rm = T)
  
  df <- data.frame(stateroute = rte, global_extract = test_rte$mean_tmax, 
                   local_extract = local_extract, local_mean = local_mean)
  
  results <- rbind(results, df)
  
}

na_points <- na_daymet %>%
  st_centroid()

daymet_points <- raster::extract(daymet_mean, na_points)

na_points$tmax_points <- daymet_points

results_points <- data.frame(stateroute = c(), global_extract = c(), local_extract = c())

for(rte in na_points$rteno) {
  test_rte <- filter(na_points, rteno == rte)
  
  local_extract <- raster::extract(daymet_mean, test_rte)
  
  df <- data.frame(stateroute = rte, global_extract = test_rte$tmax_points, 
                   local_extract = local_extract)
  
  results_points <- rbind(results_points, df)
  
}

results_polygons <- results %>%
  group_by(stateroute) %>%
  summarize(global_extract_polygon = mean(global_extract),
            local_extract_polygon = mean(local_extract),
            local_mean_polygon = mean(local_mean))

results_join <- results_points %>%
  group_by(stateroute) %>%
  summarize(global_extract_points = mean(global_extract),
            local_extract_points = mean(local_extract)) %>%
  left_join(results_polygons)

cor(results_join[, 2:6], use = "pairwise.complete.obs")
# Within Daymet: high correlation between local extracted polygons, local polygon means, global extracted points, local extracted points
# Negative correlation with global extracted polygons and the rest

# Plot Global extracted points on daymet raster

na_daymet_local <- na_daymet %>%
  left_join(results_polygons, by = c("rteno" = "stateroute"))

library(tmap)
daymet_plot <- tm_shape(daymet_mean) + tm_raster() + tm_shape(na_daymet_local) + tm_fill(col = "local_extract_polygon")
tmap_save(daymet_plot, "daymet_1990_tmax_plot.pdf")

# Do same analysis for PRISM

require(raster)
require(tidyverse)
library(prism)
library(rgdal)
library(rgeos)
library(sf)

options(prism.path = "C:/Users/gdicecco/Desktop/data/prism_1990_2018/")

setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")

us_routes <- read_sf("bbsroutes_5km_buffer.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry) %>%
  mutate(country = "US")

na_routes <- us_routes %>%
  st_transform("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

prism <- prism_stack(filter(prism_files, year == 2000, grepl("tmax", env))$files)
prism_mean <- mean(prism, na.rm = T)
routePRISM <- raster::extract(prism_mean, na_routes, fun = mean, na.rm = T, df = T) 

na_prism <- bind_cols(na_routes, routePRISM)

us_points <- na_routes %>%
  st_centroid()

routePRISM_points <- raster::extract(prism_mean, us_points, df = T)

na_prism_points <- bind_cols(us_points, routePRISM_points) %>%
  dplyr::select(rteno, layer) %>%
  st_set_geometry(NULL) %>%
  rename(prism_points = "layer")

tm_shape(prism_mean) + tm_raster() + tm_shape(na_prism) + tm_polygons(col = "layer")

# Compare these groups between PRISM and Daymet

results_compare <- results_join %>%
  left_join(na_prism, by = c("stateroute" = "rteno")) %>%
  left_join(na_prism_points, by = c("stateroute" = "rteno"))

cor(results_compare[, c(2:6, 12:13)], use = "pairwise.complete.obs")
