## Analysis of route-level environmental change (climate, habitat fragmentation)

library(tidyverse)
library(raster)
library(rgdal)
library(tmap)
library(sf)
library(ggplot2)
library(cowplot)

######## Reading in and subsetting data ##########
# Population data
## BBS
## List of species observed in BCRs of interest during time window (1990-present)
routes <- read.csv("\\\\BioArk\\HurlbertLab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

routes <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_routes_20170712.csv")
weather <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_weather_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)
setwd("\\\\BioArk/hurlbertlab/Databases/BBS/GPS_stoplocations/")
setwd("/Volumes/hurlbertlab/Databases/BBS/GPS_stoplocations/")
us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data$rteno > 3000), ]

# plot of routes

setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
setwd("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us.proj <- readOGR("BCRs_contiguous_us.shp")

us_subs_transf <- spTransform(us_subs, crs(us.proj))

plot(us.proj, col = "gray73", border = "gray73")
plot(us_subs_transf, add = T)

routes.short <- RT1.routes %>% # subset stateroutes that were filtered by criteria above
  filter(stateroute %in% us_subs@data$rteno)
# 2161 routes

# Habitat fragmentation data
frags <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F)
frags <- read.csv("/Volumes/hurlbertlab/DiCecco/data/fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F)

route_ed <- frags %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge)/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2011` - `1992`)

# Climate data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/climate/")
climate_trends <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)

climate_wide <- climate_trends %>%
  dplyr::select(-trendPval) %>%
  spread(key = "env", value = "climateTrend")

route_trends <- climate_wide %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(route_ed) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude))

####### Route level trends in environmental change #########

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/")
ggplot(route_trends, aes(x = tmax, y = deltaED)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmax", y = "Trend in deltaED")
ggsave("route_tmax_ded.tiff", units = "in")

ggplot(route_trends, aes(x = tmin, y = deltaED)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmin", y = "Trend in deltaED")
ggsave("route_tmin_ded.tiff", units = "in")

ggplot(route_trends, aes(x = tmin, y = ppt)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmin", y = "Trend in ppt")
ggsave("route_tmin_ppt.tiff", units = "in")

ggplot(route_trends, aes(x = tmax, y = ppt)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmax", y = "Trend in ppt")
ggsave("route_tmax_ppt.tiff", units = "in")

ggplot(route_trends, aes(x = tmax, y = tmin)) + geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0, cex = 1, color = "red", lty = 2) +
  geom_vline(xintercept = 0, cex = 1, color = "blue", lty = 2) +
  labs(x = "Trend in Tmax", y = "Trend in Tmin")
ggsave("route_tmax_tmin.tiff", units = "in")


######### Univariate route maps #######
# Map of route changes - univariate
## Climate trends are not z scores (see climate/script_prism_climate_trend.R)

setwd("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us_sf <- read_sf("BCRs_contiguous_us.shp")

routes_sf <- st_as_sf(route_trends, coords = c("longitude", "latitude"))

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/")
us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "gray")

tmax_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "tmax", palette = "-RdBu", midpoint = NA, size = 0.2, style = "cont", title = "Tmax")
tmax_map
tmap_save(tmax_map, "routes_tmax_map.tiff", units = "in")

tmin_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "tmin", palette = "-RdBu", midpoint = NA, size = 0.2, style = "cont", title = "Tmin")
tmin_map
tmap_save(tmin_map, "routes_tmin_map.tiff", units = "in")

ppt_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "ppt", palette = "PRGn", midpoint = NA, size = 0.2, style = "cont", title = "Ppt")
ppt_map
tmap_save(ppt_map, "routes_ppt_map.tiff", units = "in")

ded_map <- us + tm_shape(routes_sf) + 
  tm_dots(col = "deltaED", palette = "-PiYG", midpoint = NA, size = 0.2, style = "cont", title = "deltaED")
ded_map
tmap_save(ded_map, "routes_dED_map.tiff", units = "in")

######## Bivariate plots ######

# Landcover legend
newcode <- data.frame(code = seq(1,9), 
                      legend = c("Open water", "Urban", "Barren", "Forest", "Shrubland", 
                                 "Agricultural", "Grasslands", "Wetlands", "Perennial ice, snow"))

# Bivariate color labels
d<-expand.grid(x=1:3,y=1:3)
d<-merge(d,data.frame(x=1:3,xlabel=c("Δ Prop.landscape low", "Δ Prop.landscape middle","Δ Prop.landscape high")),by="x")
d<-merge(d,data.frame(y=1:3,ylabel=c("Δ Edge density low", "Δ Edge density middle","Δ Edge density high")),by="y")
d$xy <- paste0(d$x, d$y)
hex <- c("#EBF4F3", "#DFF2C4", "#FEF286", "#C2BFD4", "#8BC2BD", "#84CC8C", "#9675A0", "#6880A6", "#2C8F8A")
d$hex <- hex   

# Legend for bivariate map
g.legend<-
  ggplot(d, aes(x,y,fill=atan(y/x),alpha=x+y))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_void()+
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.margin=margin(t=10,b=10,l=10))+
  labs(x="Δ Prop.landscape",y="Δ Edge density")+
  theme(axis.title=element_text(color="black", size = 12), 
        axis.title.y = element_text(angle = 90))+
  # Draw some arrows:
  geom_segment(aes(x=1, xend = 3 , y=0, yend = 0), size=1.5,
               arrow = arrow(length = unit(0.6,"cm"))) +
  geom_segment(aes(x=0, xend = 0 , y=1, yend = 3), size=1.5,
               arrow = arrow(length = unit(0.6,"cm"))) 

#### Change in forest: fragmentation and loss

forest_ed <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  dplyr::select(stateroute, year, edge.density) %>%
  spread(key = "year", value = "edge.density") %>%
  mutate(deltaED = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))

forest_map <- us + tm_shape(forest_ed) + 
  tm_dots(col = "deltaED", palette = "-RdYlGn", midpoint = NA, size = 0.2, style = "cont", title = "deltaED Forest")
forest_map 
tmap_save(forest_map, "routes_dED_forest_map.tiff", units = "in")

forest_prop <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  dplyr::select(stateroute, year, prop.landscape) %>%
  spread(key = "year", value = "prop.landscape") %>%
  mutate(deltaPL = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))

forest_change <- forest_ed %>%
  dplyr::select(stateroute, deltaED, geometry) %>%
  st_join(dplyr::select(forest_prop, stateroute, deltaPL, geometry))

ggplot(forest_change, aes(x = deltaPL, y = deltaED)) + geom_point(alpha = 0.75) + 
  geom_hline(yintercept = 0, col = "red", lty = 2, cex = 1) +
  geom_vline(xintercept = 0, col = "blue", lty = 2, cex = 1) +
  labs(x = "Change in proportion of landscape", y = "Change in edge density", title = "Forest cover")
ggsave("forest_cover.pdf", units = "in")

x.q <- quantile(forest_change$deltaPL, c(0.33, 0.66, 1), na.rm = T)
y.q <- quantile(forest_change$deltaED, c(0.33, 0.66, 1), na.rm = T)

forest_change_biv <- forest_change %>%
  mutate(y = ifelse(deltaED < y.q[1], 1, ifelse(deltaED < y.q[2], 2, 3)), 
         x = ifelse(deltaPL < x.q[1], 1, ifelse(deltaPL < x.q[2], 2, 3))) %>%
  mutate(xy = paste0(x, y)) %>%
  left_join(dplyr::select(d, xy, hex))

ggplot(forest_change_biv, aes(x = deltaPL, y = deltaED, color = atan(y/x), alpha = x + y)) +
  geom_point(size = 1) + guides(alpha = F, color = F) +
  geom_hline(yintercept = y.q, color = "gray20", lty = 2) +
  geom_vline(xintercept = x.q, color = "gray20", lty = 2) +
  scale_color_viridis_c(name = "Color path") +
  labs(x = "Change in proportion of landscape", y = "Change in edge density", title = "Forest cover")
ggsave("forest_cover_bivariate.pdf", units = "in")

# Bivariate forest fragmentation and loss map
us <- tm_shape(us_sf) + tm_borders(col = "black") + tm_fill(col = "gray40")
forest_map <- us + tm_shape(forest_change_biv) + 
  tm_dots(col = "hex", size = 0.2)
forest_map 
print(g.legend, vp = viewport(0.22, 0.15, width = 0.25, height = 0.25))

## Change in urbanization: fragmentation and gain
forest_ed <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  dplyr::select(stateroute, year, edge.density) %>%
  spread(key = "year", value = "edge.density") %>%
  mutate(deltaED = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))

forest_map <- us + tm_shape(forest_ed) + 
  tm_dots(col = "deltaED", palette = "-RdYlGn", midpoint = NA, size = 0.2, style = "cont", title = "deltaED Forest")
forest_map 
tmap_save(forest_map, "routes_dED_forest_map.tiff", units = "in")

forest_prop <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  dplyr::select(stateroute, year, prop.landscape) %>%
  spread(key = "year", value = "prop.landscape") %>%
  mutate(deltaPL = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))

forest_change <- forest_ed %>%
  dplyr::select(stateroute, deltaED, geometry) %>%
  st_join(dplyr::select(forest_prop, stateroute, deltaPL, geometry))

ggplot(forest_change, aes(x = deltaPL, y = deltaED)) + geom_point(alpha = 0.75) + 
  geom_hline(yintercept = 0, col = "red", lty = 2, cex = 1) +
  geom_vline(xintercept = 0, col = "blue", lty = 2, cex = 1) +
  labs(x = "Change in proportion of landscape", y = "Change in edge density", title = "Forest cover")
ggsave("forest_cover.pdf", units = "in")

x.q <- quantile(forest_change$deltaPL, c(0.33, 0.66, 1), na.rm = T)
y.q <- quantile(forest_change$deltaED, c(0.33, 0.66, 1), na.rm = T)

forest_change_biv <- forest_change %>%
  mutate(y = ifelse(deltaED < y.q[1], 1, ifelse(deltaED < y.q[2], 2, 3)), 
         x = ifelse(deltaPL < x.q[1], 1, ifelse(deltaPL < x.q[2], 2, 3))) %>%
  mutate(xy = paste0(x, y)) %>%
  left_join(dplyr::select(d, xy, hex))

ggplot(forest_change_biv, aes(x = deltaPL, y = deltaED, color = atan(y/x), alpha = x + y)) +
  geom_point(size = 1) + guides(alpha = F, color = F) +
  geom_hline(yintercept = y.q, color = "gray20", lty = 2) +
  geom_vline(xintercept = x.q, color = "gray20", lty = 2) +
  scale_color_viridis_c(name = "Color path") +
  labs(x = "Change in proportion of landscape", y = "Change in edge density", title = "Forest cover")
ggsave("forest_cover_bivariate.pdf", units = "in")

# Bivariate forest fragmentation and loss map
us <- tm_shape(us_sf) + tm_borders(col = "black") + tm_fill(col = "gray40")
forest_map <- us + tm_shape(forest_change_biv) + 
  tm_dots(col = "hex", size = 0.2)
forest_map 
print(g.legend, vp = viewport(0.22, 0.15, width = 0.25, height = 0.25))

###### Grassland change


##### Urban change
urban_ed <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Urban") %>%
  dplyr::select(stateroute, year, edge.density) %>%
  spread(key = "year", value = "edge.density") %>%
  mutate(deltaED = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))

urban_prop <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Urban") %>%
  dplyr::select(stateroute, year, prop.landscape) %>%
  spread(key = "year", value = "prop.landscape") %>%
  mutate(deltaPL = `2011` - `1992`) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"))

urban_change <- urban_ed %>%
  dplyr::select(stateroute, deltaED, geometry) %>%
  st_join(dplyr::select(urban_prop, stateroute, deltaPL, geometry))

ggplot(urban_change, aes(x = deltaPL, y = deltaED)) + geom_point(alpha = 0.75) + 
  geom_hline(yintercept = 0, col = "red", lty = 2, cex = 1) +
  geom_vline(xintercept = 0, col = "blue", lty = 2, cex = 1) +
  labs(x = "Change in proportion of landscape", y = "Change in edge density", title = "Urban cover")
ggsave("urban_cover.pdf", units = "in")

x.q <- quantile(urban_change$deltaPL, c(0.33, 0.66, 1), na.rm = T)
y.q <- quantile(urban_change$deltaED, c(0.33, 0.66, 1), na.rm = T)

urban_change_biv <- urban_change %>%
  mutate(y = ifelse(deltaED < y.q[1], 1, ifelse(deltaED < y.q[2], 2, 3)), 
         x = ifelse(deltaPL < x.q[1], 1, ifelse(deltaPL < x.q[2], 2, 3))) %>%
  mutate(xy = paste0(x, y)) %>%
  left_join(dplyr::select(d, xy, hex))

ggplot(urban_change_biv, aes(x = deltaPL, y = deltaED, color = atan(y/x), alpha = x + y)) +
  geom_point(size = 1) + guides(alpha = F, color = F) +
  geom_hline(yintercept = y.q, color = "gray20", lty = 2) +
  geom_vline(xintercept = x.q, color = "gray20", lty = 2) +
  scale_color_viridis_c(name = "Color path") +
  labs(x = "Change in proportion of landscape", y = "Change in edge density", title = "Urban cover")
ggsave("urban_cover_bivariate.pdf", units = "in")

# Bivariate urban fragmentation and loss map
us <- tm_shape(us_sf) + tm_borders(col = "black") + tm_fill(col = "gray40")
urban_map <- us + tm_shape(na.omit(urban_change_biv)) + 
  tm_dots(col = "hex", size = 0.2)
urban_map 
print(g.legend, vp = viewport(0.22, 0.15, width = 0.25, height = 0.25))

#### Land cover plus climate maps

##### Moran's I for env variables ########

##### Map of route-level abundance trends #######
# Subset species: diurnal land birds
landbirds <- species %>%
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010) %>%
  filter(aou != 22860) # Eurasian collared dove

## Population trends
counts.subs <- counts %>%
  filter(aou %in% landbirds$aou) %>%
  merge(routes.short, by = c("stateroute", "year")) %>%
  filter(year > 1990, year < 2017)
# 2031 routes

library(purrr)
library(broom)
abund_trend <- counts.subs %>%
  group_by(stateroute) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    df.short <- df %>%
      group_by(year) %>%
      summarize(abund = sum(speciestotal)) %>%
      dplyr::select(year, abund) %>%
      unique()
    lm(abund ~ year, df.short)
  })) %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    length(unique(df$year))
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  mutate(abundTrend = map_dbl(lm_broom, ~{
    df <- .
    df$estimate[2]
  })) %>%
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    df$p.value[2]
  })) 

abund_sf <- abund_trend %>%
  filter(nObs > 9) %>%
  dplyr::select(-data, -lmFit, -lm_broom) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>% 
  filter(abundTrend < 200) %>% # a couple outliers in LA
  st_as_sf(coords = c("longitude", "latitude"))

bird_map <- us + tm_shape(abund_sf) + 
  tm_dots(col = "abundTrend", palette = "-RdBu", size = 0.2, style = "cont", title = "Abundance trend")
bird_map
tmap_save(bird_map, "routes_birdAbund_map.tiff", units = "in")

              