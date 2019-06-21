## Make figures for cartoon methods

# Rasters of tmin, tmax, ppt for one route
# Example routes with ~0 ED, 0.25 ED, 0.5 ED, 1 ED

require(raster)
require(tidyverse)
library(prism)
library(rgdal)
library(rgeos)
library(tmap)
library(cowplot)

options(prism.path = "/Users/gracedicecco/Desktop/prism_2018/")

# Get the PRISM data (1 time only)
# Annual precip, monthly temperature mins and maxes for breeding season
get_prism_annual("ppt", years = 2017, keepZip = F)

get_prism_monthlys("tmin", years = 2017, mon = c(5), keepZip = F)

get_prism_monthlys("tmax", years = 2017, mon = c(5), keepZip = F)

# Extract breeding season average temps 
prism_files <- ls_prism_data()
prism_files$env <- word(prism_files$files, start = 2, sep = "_")
prism_files$date <- word(prism_files$files, start = 5, sep = "_")
prism_files$year <- substr(prism_files$date, 1, 4)
prism_files$month <- substr(prism_files$date, 5, 6)

prism <- prism_stack(filter(prism_files, year != 2018)$files) # leave out provisional 2018 data for now
prismCRS <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")
bufferRoutes <- readOGR("bbsroutes_5km_buffer.shp")
bufferRoutes_transf <- spTransform(bufferRoutes, prismCRS)

route1 <- subset(bufferRoutes_transf, rteno == 2001)

list2env(setNames(unstack(prism), names(prism)), .GlobalEnv)

setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/figures/")
prism_crop <- crop(PRISM_ppt_stable_4kmM3_2017_bil, route1)
route_prism <- mask(prism_crop, route1)
tm_shape(route_prism) +
  tm_raster(style = "cont", title = "Ppt", palette = "YlGn")
ggsave("methods_ppt_raster.tiff", units = "in")

prism_crop_tmax <- crop(PRISM_tmax_stable_4kmM2_201705_bil, route1)
route_prism_tmax <- mask(prism_crop_tmax, route1)
tm_shape(route_prism_tmax) +
  tm_raster(style = "cont", title = "Tmax", palette = "Reds", legend.show = T)

prism_crop_tmin <- crop(PRISM_tmin_stable_4kmM2_201705_bil, route1)
route_prism_tmin <- mask(prism_crop_tmin, route1)
tm_shape(route_prism_tmax) +
  tm_raster(style = "cont", title = "Tmin", palette = "Blues")

# Plots of abundance trend and climate trend for one species

setwd("/Volumes/hurlbertlab/DiCecco/data/")
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
routePRISM <- read.csv("bbs_routes_breeding_season_climate.csv", stringsAsFactors = F)
routeClim <- filter(routePRISM, ID == 6) %>%
  group_by(year, env) %>%
  summarize(mean = mean(val))

library(ggplot2)
theme_set(theme_classic())

ggplot(filter(routeClim, env == "ppt"), aes(x = year, y = mean)) + geom_point(size = 4) + 
  geom_smooth(method = "lm", se = F, color = "green", cex = 3) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.line = element_line(colour = 'black', size = 3)) +
  labs(x = "Year", y = "Precipitation (mm)")

ggplot(filter(routeClim, env == "tmax"), aes(x = year, y = mean)) + geom_point(size = 4) + 
  geom_smooth(method = "lm", se = F, color = "red", cex = 3) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.line = element_line(colour = 'black', size = 3)) +
  labs(x = "Year", y = "Tmax (C)")

ggplot(filter(routeClim, env == "tmin"), aes(x = year, y = mean)) + geom_point(size = 4) + 
  geom_smooth(method = "lm", se = F, cex = 3) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.line = element_line(colour = 'black', size = 3)) +
  labs(x = "Year", y = "Tmin (C)")

routes <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_routes_20170712.csv")
counts <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_counts_20170712.csv")
species <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_species_20170712.csv")
weather <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_weather_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

counts_onespp <- counts %>%
  filter(aou == 2890, year > 1989, stateroute == 2001)

ggplot(counts_onespp, aes(x = year, y = speciestotal)) + geom_point(size = 4) + 
  geom_smooth(method = "lm", se = F, color = "purple", cex = 3) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.line = element_line(colour = 'black', size = 3)) +
  labs(x = "Year", y = "Abundance")

### Plot of example routes for edge density scale
setwd("//BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")
setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/")
bufferRoutes <- readOGR("bbsroutes_5km_buffer.shp")

nlcd2016 <- raster("C:/Users/gdicecco/Desktop/data/nlcd/NLCD_2016_Land_Cover_L48_20190424/nlcd_2016_whole_simplified.tif")

frags <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F)
frags <- read.csv("/Volumes/hurlbertlab/DiCecco/data/fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F)

newcode <- data.frame(code = seq(1,9), 
                      legend = c("Open water", "Urban", "Barren", "Forest", "Shrubland", 
                                 "Agricultural", "Grasslands", "Wetlands", "Perennial ice, snow"))

## landscape edge density

route_ed <- frags %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge)/sum(total.area)) %>%
  filter(year == 2016) %>% # use filter to ID routes
  arrange(ED)

bufferRoutes_transf <- spTransform(bufferRoutes, crs(nlcd2011))

route0.01 <- subset(bufferRoutes_transf, rteno == 55003) #0.014
nlcd_crop <- crop(nlcd2011, route0.01)
ED0.01 <- mask(nlcd_crop, route0.01)

route0.251 <- subset(bufferRoutes_transf, rteno == 34305)
nlcd_crop0.25 <- crop(nlcd2011, route0.251)
ED0.251 <- mask(nlcd_crop0.25, route0.251)

route0.5 <- subset(bufferRoutes_transf, rteno == 46052)
nlcd_crop0.5 <- crop(nlcd2011, route0.5)
ED0.5 <- mask(nlcd_crop0.5, route0.5)

route1 <- subset(bufferRoutes_transf, rteno == 92039)
nlcd_crop1 <- crop(nlcd2011, route1)
ED1.0 <- mask(nlcd_crop1, route1)

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
pdf("figures/methods_figs/edge_density_scale.pdf", height = 6, width = 8)
par(mfrow = c(2,2))
plot(ED0.01, main = "Edge density = 0.014")
plot(ED0.251, main = "Edge density = 0.251")
plot(ED0.5, main = "Edge density = 0.500")
plot(ED1.0, main = "Edge density = 1.000")
dev.off()

## forest edge density

forest_ed <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  mutate(sum.area = sum(total.area)) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(ED = total.edge/sum.area,
            propForest = prop.landscape) %>%
  filter(year == 2016, propForest> 0.2) %>% # use filter to ID routes
  arrange(ED)

bufferRoutes_transf <- spTransform(bufferRoutes, crs(nlcd2016))

routeLow <- subset(bufferRoutes_transf, rteno == 63908)
nlcd_crop <- crop(nlcd2016, routeLow)
EDlow <- mask(nlcd_crop, routeLow)
# ED = 0.059, % for = 0.965
# Busick, NC

routeMed <- subset(bufferRoutes_transf, rteno == 85106)
nlcd_cropMed <- crop(nlcd2016, routeMed)
EDmed <- mask(nlcd_cropMed, routeMed)
# ED = 0.194, % for = 0.695
# Uinta National Forest, UT

routeHigh <- subset(bufferRoutes_transf, rteno == 18010)
nlcd_cropHigh <- crop(nlcd2016, routeHigh)
EDHigh <- mask(nlcd_cropHigh, routeHigh)
# ED = 0.328, % for = 0.529
# Greenwich, CT

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/")

#### Figure 1 in manuscript - edge density scale
nlcd_palette <- c("#0000FF", "#FF9900", "#E5E5CC", "#006600", "#B2B200", "#FFB3CC", "#E5CC99", "#80FFCC", "#FFFFFF")
low <- tm_shape(EDlow) + tm_raster(palette = nlcd_palette, style = "cat", labels = as.character(newcode$legend), title = "") +
  tm_scale_bar(size = 2, breaks = c(0, 5)) + 
  tm_layout(legend.text.size = 1, title = "ED = 0.059\nForest cover = 0.965", 
            title.position = c("left", "bottom"),
            title.size = 1,
            main.title = "A",
            main.title.fontface = "bold")
tmap_save(low, "figures/methods_figs/forest_ed_low.tiff")

med <- tm_shape(EDmed) + tm_raster(palette = nlcd_palette, style = "cat") +
  tm_legend(show = F) + tm_scale_bar(size = 2, breaks = c(0, 5)) + 
  tm_layout(legend.text.size = 1, title = "ED = 0.194\nForest cover = 0.695", 
            title.position = c("left", "bottom"),
            title.size = 1,
            inner.margins = c(0, 0.13, 0, 0.13),
            main.title = "B",
            main.title.fontface = "bold")
tmap_save(med, "figures/methods_figs/forest_ed_med.tiff")

high <- tm_shape(EDHigh) + tm_raster(palette = nlcd_palette, style = "cat") +
  tm_legend(show = F)+ tm_scale_bar(size = 2, position = c("RIGHT", "BOTTOM"), breaks = c(0, 5)) + 
  tm_layout(legend.text.size = 1, title = "ED = 0.328\nForest cover = 0.529", 
            title.size = 1,
            title.position = c("right", "bottom"),             
            main.title = "C",
            main.title.fontface = "bold")
tmap_save(high, "figures/methods_figs/forest_ed_high.tiff")

arrange <- tmap_arrange(low, med, high, ncol = 3)
tmap_save(arrange, "figures/methods_figs/forest_ed_multipanel.pdf", units = "in", height = 6, width = 12)
