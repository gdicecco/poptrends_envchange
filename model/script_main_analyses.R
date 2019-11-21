### Fit species-specific models using area-sensitivity classification

#### Libraries ####

library(tidyverse)
library(raster)
library(rgdal)
library(purrr)
library(broom)
library(tmap)
library(cowplot)
library(sf)
library(forcats)
library(grid)
library(spdep)
library(ncf)

### Read in data #####

# Bird population data
## BBS
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
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# species four letter codes
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
spp_codes <- read.csv("four_letter_codes_birdspp.csv", stringsAsFactors = F)
spp_codes$common_name_lower <- tolower(spp_codes$COMMONNAME)

# subset routes that are between 38000 and 42000 m, no Alaska
setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/")

ca_routes <- read_sf("DiCecco/nlcd_frag_proj_shapefiles/bbs_canada_route_areas.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry)
us_routes <- read_sf("DiCecco/nlcd_frag_proj_shapefiles/bbsroutes_5km_buffer.shp") %>%
  dplyr::select(rteno, RTENAME, STATUS, geometry)

na_routes <- rbind(ca_routes, us_routes)

na_route_paths <- read_sf("Databases/BBS/GPS_stoplocations/BBS_routes_USandCanada/bbs_routes_usandcanada.shp")

# Subset stateroutes: runtype = 1, sampled once every 5 years from 1992-2016 in US and 1990-2010 in Canada, 
# between 38 and 42 km in length

routes.us <- RT1.routes %>%
  filter(stateroute %in% na_routes$rteno) %>%
  mutate(year_bin = 5*floor(year/5) + 5/2) %>%
  filter(countrynum == 840) %>%
  filter(year >= 1992 & year <= 2016) %>%
  group_by(stateroute) %>%
  distinct(year_bin) %>%
  count() %>%
  filter(n >= 5)

routes.ca <- RT1.routes %>%
  filter(stateroute %in% na_routes$rteno) %>%
  mutate(year_bin = 5*floor(year/5) + 5/2) %>%
  filter(countrynum == 124) %>%
  filter(year >= 1990 & year <= 2010) %>%
  group_by(stateroute) %>%
  distinct(year_bin) %>%
  count() %>%
  filter(n >= 4)

routes.short <- RT1.routes %>%
  filter(stateroute %in% routes.us$stateroute | stateroute %in% routes.ca$stateroute)
# 1834 routes

# Subset species: forest breeding birds

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation")
area_sensitive <- read.csv("traits/area-sensitivity-forest-birds-Bouliner1998.csv", stringsAsFactors = F)
area_sensitive <- area_sensitive[-12, -c(4)]

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
  left_join(spp_codes)

fourletter_codes <- dplyr::select(area_aous, aou, SPEC)

## Population trends
counts.subs <- counts %>%
  filter(aou %in% area_aous$aou) %>%
  merge(filter(RT1.routes, stateroute %in% routes.short$stateroute), by = c("stateroute", "year")) %>%
  filter(rpid == 101) %>%
  filter(year >= 1990, year < 2017)
# 1834 routes

# # first year observers
# 
# first_year <- function(obsns) {
#   first_yrs <- c()
#   for(i in 1:length(obsns)) {
#     ob <- obsns[i]
#     obs <- obsns[1:i]
#     first_yrs <- c(first_yrs, ifelse(length(obs[obs == ob]) == 1, 1, 0))
#   }
#   return(first_yrs)
# }
# 
# obs_years <- weather %>%
#   group_by(stateroute) %>%
#   arrange(year) %>%
#   nest() %>%
#   mutate(first_yr = purrr::map(data, ~{
#     df <- .
#     first_year(df$obsn)
#   })) %>%
#   unnest() %>%
#   dplyr::select(stateroute, year, first_yr)
# 
# # abundance trends read in
# # 1990-2010 for Canada, 1992-2016 for US
# abund_trend <- counts.subs %>%
#   left_join(obs_years, by = c('stateroute', 'year')) %>%
#   group_by(aou, stateroute) %>%
#   nest() %>%
#   mutate(nObs = map_dbl(data, ~{
#     df <- .
#     country <- unique(df$countrynum.x)
#     if(country == 840){
#       df <- df %>%
#         filter(year >= 1992 & year <= 2016) %>%
#         unique()
#       nrow(df)
#     } else {
#       df <- df %>%
#         filter(year >= 1990 & year <= 2010) %>%
#         unique()
#       nrow(df)
#     }
#   })) %>%
#   filter(nObs > 9) %>%
#   mutate(lmFit = purrr::map(data, ~{
#     df <- .
#     country <- unique(df$countrynum.x)
#     if(country == 840){
#       df.short <- df %>%
#         dplyr::select(year, first_yr, speciestotal) %>%
#         filter(year >= 1992 & year <= 2016) %>%
#         unique()
#       glm(speciestotal ~ year + first_yr, family = poisson, data = df.short)
#     } else {
#       df.short <- df %>%
#         dplyr::select(year, first_yr, speciestotal) %>%
#         filter(year >= 1990 & year <= 2010) %>%
#         unique()
#       glm(speciestotal ~ year + first_yr, family = poisson, data = df.short)
#     }
#   }))  %>%
#   mutate(lm_broom = map(lmFit, tidy)) %>%
#   mutate(abundTrend = map_dbl(lm_broom, ~{
#     df <- .
#     df$estimate[2]
#   })) %>%
#   mutate(trendInt = map_dbl(lm_broom, ~{
#     df <- .
#     df$estimate[1]
#   })) %>%
#   mutate(trendPval = map_dbl(lm_broom, ~{
#     df <- .
#     df$p.value[2]
#   })) %>%
#   mutate(obsTrend = map_dbl(lm_broom, ~{
#     df <- .
#     df$estimate[3]
#   })) %>%
#   mutate(obsPval = map_dbl(lm_broom, ~{
#     df <- .
#     df$p.value[3]
#   }))
# 
# hist(abund_trend$abundTrend)
# hist(abund_trend$nObs)
# 
# # Save abundance trend data
# abund_trend <- abund_trend %>%
#   dplyr::select(-data, -lmFit, -lm_broom)
# 
# setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
# write.csv(abund_trend, "BBS_abundance_trends.csv", row.names = F)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
abund_trend <- read.csv("BBS_abundance_trends.csv", stringsAsFactors = F) %>%
  filter(aou %in% area_aous$aou) %>% ## Only ones lost from area_aous are nocturnal: owls, nightjars (Caprimulgus)
  left_join(area_aous) %>%
  filter(sporder != "Accipitriformes" & sporder != "Strigiformes" & sporder != "Caprimulgiformes") ## Remove hawks

length(abund_trend$trendPval[abund_trend$trendPval < 0.05])/nrow(abund_trend) # 45%
length(abund_trend$trendPval[abund_trend$trendPval < 0.01])/nrow(abund_trend) # 34%
length(abund_trend$trendPval[abund_trend$trendPval < 0.0001])/nrow(abund_trend) # 20%

abund_trend_sig <- abund_trend %>%
  mutate(abund_class = case_when(trendPval >= 0.05 ~ "p > 0.05",
                                 trendPval >= 0.01 & trendPval < 0.05 ~ "0.01 < p < 0.05",
                                 trendPval >= 0.001 & trendPval < 0.01 ~ "0.001 < p < 0.01",
                                 trendPval < 0.001 ~ "p < 0.001"),
         obs_class = case_when(obsPval >= 0.05 ~ "p > 0.05",
                               obsPval >= 0.01 & obsPval < 0.05 ~ "0.01 < p < 0.05",
                               obsPval >= 0.001 & obsPval < 0.01 ~ "0.001 < p < 0.01",
                               obsPval < 0.001 ~ "p < 0.001"))

# Trait data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/")
traits <- read.csv("traits/spp_traits.csv", stringsAsFactors = F)

traits.short <- traits %>%
  dplyr::select(Common_name, aou, nHabitats1, nHabitats2, volume)

setwd("//BioArk/HurlbertLab/DiCecco/Data/")
setwd("/Volumes/hurlbertlab/dicecco/data/")
correlates <- read.csv("Master_RO_Correlates_20110610.csv", stringsAsFactors = F)

breedingRange <- correlates %>%
  dplyr::select(AOU, Brange_Area_km2)

# Climate data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/climate/")
climate_trends <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)

climate_wide <- climate_trends %>%
  dplyr::select(-trendPval) %>%
  spread(key = "env", value = "climateTrend")

# Habitat fragmentation data
setwd("/Volumes/hurlbertlab/DiCecco/data/")
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
frags <- read.csv("fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F)
frags_ca <- read.csv("fragmentation_indices_canada.csv", stringsAsFactors = F) %>%
  group_by(year, stateroute) %>% 
  mutate(n_zones = as.factor(n_distinct(file))) %>%
  filter(n_zones == 1)

# Landcover legend US
newcode <- data.frame(code = seq(1,9), 
                      legend = c("Open water", "Urban", "Barren", "Forest", "Shrubland", 
                                 "Agricultural", "Grasslands", "Wetlands", "Perennial ice, snow"))

# Landcover legend Canada
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
canada_code <- read.csv("landcover-finescale/scripts_canada_landcover/canada_landcover_classification.csv", stringsAsFactors = F)
colnames(canada_code)[1:2] <- c("legend", "code")

# Forest fragmentation deltas for US
# long to wide

forest_ed <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  mutate(sum.area = sum(total.area)) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(ED = total.edge/sum.area,
            propForest = prop.landscape,
            meanPatchArea = mean.patch.area) %>%
  filter(year == 2016, propForest > 0.25) %>% # use filter to ID routes
  arrange(ED)

forest_deltaED <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge[legend == "Forest"])/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2016` - `1992`,
            deltaEDrecent = `2016` - `2013`,
            recentChange = abs(deltaEDrecent)/abs(deltaED)) %>%
  filter(recentChange < 0.5)

forest_deltaP <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(propForest = prop.landscape) %>%
  spread(key = "year", value = "propForest") %>%
  group_by(stateroute) %>%
  summarize(deltaProp = `2016` - `1992`)

# Forest fragmentation deltas for Canada
# long to wide

forest_ed_ca <- frags_ca %>%
  left_join(canada_code, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  mutate(sum.area = sum(total.area)) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(ED = total.edge/sum.area,
            propForest = prop.landscape,
            meanPatchArea = mean.patch.area) %>%
  filter(year == 2010, propForest > 0.25) %>% # use filter to ID routes
  arrange(ED)

forest_deltaED_ca <- frags_ca %>% 
  left_join(canada_code, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge[legend == "Forest"])/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2010` - `1990`,
            deltaEDrecent = `2010` - `2000`,
            recentChange = abs(deltaEDrecent)/abs(deltaED)) %>%
  filter(recentChange < 0.5)

forest_deltaP_ca <- frags_ca %>%
  left_join(canada_code, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(propForest = prop.landscape) %>%
  spread(key = "year", value = "propForest") %>%
  group_by(stateroute) %>%
  summarize(deltaProp = `2010` - `1990`)


# join edge density, forest cover for US and canada

forest_us <- forest_ed %>%
  right_join(forest_deltaED) %>%
  left_join(forest_deltaP) %>%
  mutate_at(c("deltaED", "deltaProp"), function(x) x/25) # divide change in ED and proportion cover by number of years in dataset

forest_ca <- forest_ed_ca %>%
  right_join(forest_deltaED_ca) %>%
  left_join(forest_deltaP_ca) %>%
  mutate_at(c("deltaED", "deltaProp"), function(x) x/21) # divide change in ED and proportion cover by number of years in dataset

forest <- bind_rows(forest_us, forest_ca)

# Master data table
env_change <- forest %>%
  left_join(climate_wide, by = "stateroute")
# write.csv(env_change, "model/bbs_route_env_change.csv", row.names = F)

clim_hab_poptrend <- abund_trend %>%
  left_join(dplyr::select(traits.short, -Common_name)) %>%
  left_join(forest, by = "stateroute") %>%
  left_join(climate_wide, by = "stateroute") %>%
  filter(!is.na(propForest))

# shapefile of US/Canada
map <- read_sf("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/landcover-finescale/ne_50m_admin_1_states_provinces_lakes.shp")
na_map <- map %>%
  filter(sr_adm0_a3 %in% c("CAN", "USA")) %>%
  st_crop(c(ymin = 20, ymax = 60, xmax = -52, xmin = -145))

##### Analysis #####
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/")

#### Route trends ####

# Plot of routes used

forest_routes <- na_route_paths %>%
  filter(rteno %in% forest_ed$stateroute | rteno %in% forest_ed_ca$stateroute)

forest_transf <- st_transform(forest_routes, crs(na_map))

hist(clim_hab_poptrend$deltaED)
hist(clim_hab_poptrend$deltaProp)

cor(clim_hab_poptrend$ED, clim_hab_poptrend$propForest) # -0.32
cor(clim_hab_poptrend$deltaED, clim_hab_poptrend$deltaProp) # -0.22
cor(dplyr::select(clim_hab_poptrend, tmin, tmax, deltaED, deltaProp, ED, propForest), use = "pairwise.complete.obs")

# Figure: map of routes

color_scale <- data.frame(color = c(1:4), 
                          dED_hex = c("#5aae61", "#e6f5d0", "#fde0ef", "#de77ae"),
                          dProp_hex = c("#7b3294", "#c2a5cf", "#d9f0d3", "#5aae61"),
                          temp_hex = c("#92C5DE", "#FDDBC7", "#EF8A62", "#B2182B"), stringsAsFactors = F)

bbs_routes_forest <- forest_transf %>%
  st_as_sf() %>%
  left_join(forest, by = c("rteno" = "stateroute")) %>%
  left_join(climate_wide, by = c("rteno" = "stateroute")) %>%
  filter(!is.na(deltaProp))

bbs_routes_forcov <- tm_shape(na_map) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + 
  tm_lines(col = "deltaProp", scale = 3, breaks = quantile(bbs_routes_forest$deltaProp), palette = color_scale$dProp_hex) +
  tm_layout(title = "A", title.fontface = "bold", legend.show = F)

bbs_routes_fored <- tm_shape(na_map) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + 
  tm_lines(col = "deltaED", scale = 3, breaks = quantile(bbs_routes_forest$deltaED), palette = color_scale$dED_hex) +
  tm_layout(legend.show = F, title = "B", title.fontface = "bold")

bbs_routes_tmin <- tm_shape(na_map) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + 
  tm_lines(col = "tmin", scale = 3, breaks = quantile(bbs_routes_forest$tmin, na.rm = T), palette = color_scale$temp_hex) +
  tm_layout(legend.show = F, title = "C", title.fontface = "bold")

bbs_routes_tmax <- tm_shape(na_map) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + 
  tm_lines(col = "tmax", scale = 3, breaks = quantile(bbs_routes_forest$tmax, na.rm = T), palette = color_scale$temp_hex) +
  tm_layout(legend.show = F, title = "D", title.fontface = "bold")

bbs_routes <- tmap_arrange(bbs_routes_forcov, bbs_routes_fored, bbs_routes_tmin, bbs_routes_tmax, nrow = 2)


# Histogram legends for each variable

# Add color scale

dED_quant <- quantile(bbs_routes_forest$deltaED)
dProp_quant <- quantile(bbs_routes_forest$deltaProp)
tmin_quant <- quantile(bbs_routes_forest$tmin, na.rm = T)
tmax_quant <- quantile(bbs_routes_forest$tmax, na.rm = T)

quantile_group <- function(quantiles, value) {
  case_when(
  value <= quantiles[2] ~ 1,
  value > quantiles[2] & value <= quantiles[3] ~ 2,
  value > quantiles[3] & value <= quantiles[4] ~ 3,
  value > quantiles[4] & value <= quantiles[5] ~ 4)
}

bbs_routes_forest_plots <- bbs_routes_forest %>%
  mutate(dED_color = quantile_group(dED_quant, deltaED),
  dProp_color = quantile_group(dProp_quant, deltaProp),
  tmin_color = quantile_group(tmin_quant, tmin),
  tmax_color = quantile_group(tmax_quant, tmax))

ED_hist <- ggplot(bbs_routes_forest_plots, aes(x = deltaED, fill = as.factor(dED_color))) + geom_histogram(bins = 50) +
  labs(title = "Change in edge density") +
  scale_y_continuous(breaks = c(0, 75, 150)) + scale_x_continuous(breaks = c(-0.01, 0, 0.007)) +
  scale_fill_manual(values = color_scale$dED_hex) +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

for_hist <- ggplot(bbs_routes_forest_plots, aes(x = deltaProp, fill = as.factor(dProp_color))) + geom_histogram(bins = 30) +
  labs(title = "Change in forest cover") +
  scale_y_continuous(breaks = c(0, 100, 200)) +
  scale_fill_manual(values = color_scale$dProp_hex) +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

tmin_hist <- ggplot(bbs_routes_forest_plots, aes(x = tmin, fill = as.factor(tmin_color))) + geom_histogram(bins = 30) +
  labs(title = "Trend in Tmin") +
  scale_y_continuous(breaks = c(0, 100, 200)) +
  scale_fill_manual(values = color_scale$temp_hex) +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

tmax_hist <- ggplot(bbs_routes_forest_plots, aes(x = tmax, fill = as.factor(tmax_color))) + geom_histogram(bins = 30) +
  labs(title = "Trend in Tmax") +
  scale_y_continuous(breaks = c(0, 100, 200)) +
  scale_fill_manual(values = color_scale$temp_hex) +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10),
        axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

## Insets

vp_for <- viewport(0.1, 0.59, width = 0.175, height = 0.15)
vp_ed <- viewport(0.6, 0.59, width = 0.175, height = 0.15)

vp_tmin <- viewport(0.1, 0.09, width = 0.175, height = 0.15)
vp_tmax <- viewport(0.6, 0.09, width = 0.175, height = 0.15)

pdf("figures/methods_figs/bbs_routes.pdf", height = 10, width = 12)
print(bbs_routes)
print(for_hist, vp = vp_for)
print(ED_hist, vp = vp_ed)
print(tmin_hist, vp = vp_tmin)
print(tmax_hist, vp = vp_tmax)
dev.off()

# Appendix figure: correlation between env variables
matrix <- climate_wide %>%
  left_join(forest) %>%
  filter(!is.na(propForest)) %>%
  dplyr::select(-meanPatchArea, -stateroute, -ED, -propForest, -year, -deltaEDrecent, -recentChange) %>%
  rename("Tmin" = "tmin", "Tmax" = "tmax", "forestED" = "deltaED", "forestCover" = "deltaProp")

corr.table <- cor(matrix)
pdf("figures/methods_figs/env_cor_matrix.pdf")
corrplot::corrplot(corr.table, method = "circle", diag = F, tl.col = "black",
                   tl.cex = 1.5, cl.cex = 1)
dev.off()

# Appendix figure: Moran's I

route_trends_forest <- climate_wide %>%
  left_join(forest) %>%
  left_join(routes)

# Forest fragmentation

for.cor <- correlog(route_trends_forest$longitude, route_trends_forest$latitude, route_trends_forest$deltaED,
                    increment = 250, latlon = T, na.rm = T)

# Forest proportion cover

cov.cor <- correlog(route_trends_forest$longitude, route_trends_forest$latitude, route_trends_forest$deltaProp,
                    increment = 250, latlon = T, na.rm = T)

# Trend in tmax

tmax.cor <- correlog(route_trends_forest$longitude, route_trends_forest$latitude,route_trends_forest$tmax,
                     increment = 250, latlon = T, na.rm = T)

# Trend in tmin

tmin.cor <- correlog(route_trends_forest$longitude, route_trends_forest$latitude, route_trends_forest$tmin,
                     increment = 250, latlon = T, na.rm = T)

## Combine on one plot

env.cor <- data.frame(mean.of.class = for.cor$mean.of.class, 
                      forestCover = cov.cor$correlation,
                      forestED = for.cor$correlation,
                      tmin = tmin.cor$correlation,
                      tmax = tmax.cor$correlation) %>%
  filter(mean.of.class > 0, mean.of.class < 4100)

env.long <- tidyr::gather(env.cor[, 1:5], key = variable, value = correlation, 2:5)

ggplot(env.long, aes(x = mean.of.class, y = correlation, col = variable)) + 
  geom_point(size = 2) + geom_line(cex = 1) + geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Mean of distance class (km)", y = "Moran's I", col = "Trend") +
  scale_color_viridis_d() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
                                  legend.text = element_text(size = 12), legend.title = element_text(size = 12))
ggsave("figures/methods_figs/moransI_allenv.pdf", units = "in")

#### Species traits ####
# Breadth of forest ED and forest cover for area and non-area sensitive species

env_breadth <- clim_hab_poptrend %>%
  dplyr::group_by(Area_sensitivity, SPEC, aou) %>%
  summarize(ed = mean(ED),
            std_ed = sd(ED),
            propFor = mean(propForest),
            std_for = sd(propForest),
            patchArea = mean(meanPatchArea),
            std_patch = sd(meanPatchArea))

# Use all routes for cover_breadth

## Breadth of forest cover for all routes
#### ADD CANADIAN ROUTES
forest_ed_allroutes <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  mutate(sum.area = sum(total.area)) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(ED = total.edge/sum.area,
            propForest = prop.landscape,
            meanPatchArea = mean.patch.area) %>%
  filter(year == 2016)

forest_ed_canada <- frags_ca %>%
  left_join(canada_code, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  mutate(sum.area = sum(total.area)) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(ED = total.edge/sum.area,
            propForest = prop.landscape,
            meanPatchArea = mean.patch.area) %>%
  filter(year == 2010)

forest_allroutes <- forest_ed_allroutes %>%
  bind_rows(forest_ed_canada)

clim_hab_pop_allroutes <- abund_trend %>%
  left_join(dplyr::select(traits.short, -Common_name)) %>%
  left_join(forest_allroutes, by = "stateroute") %>%
  left_join(climate_wide, by = "stateroute") %>%
  filter(!is.na(propForest))

env_breadth_allroutes <- clim_hab_pop_allroutes %>%
  dplyr::group_by(Area_sensitivity, SPEC, aou) %>%
  summarize(propFor = mean(propForest),
            min_for = quantile(propForest, c(0.05))[[1]],
            max_for = quantile(propForest, c(0.95))[[1]],
            patchArea = mean(meanPatchArea),
            std_patch = sd(meanPatchArea))

volume <- traits %>%
  filter(aou %in% env_breadth_allroutes$aou) %>%
  left_join(dplyr::select(env_breadth_allroutes, SPEC, aou)) %>%
  dplyr::select(SPEC, aou, volume)

spp_breadths <- env_breadth_allroutes %>%
  dplyr::select(SPEC, aou, propFor, min_for, max_for) %>%
  left_join(dplyr::select(env_breadth, SPEC, aou, ed, std_ed)) %>%
  left_join(dplyr::select(traits, aou, volume))

# Species table for MS supplement

spp_table_traits <- correlates %>%
  dplyr::select(AOU, CommonName, migclass, Foraging) %>%
  left_join(spp_breadths, by = c("AOU" = "aou")) %>%
  filter(!is.na(SPEC), SPEC != "CERW") %>%
  ungroup() %>%
  dplyr::select(AOU, CommonName, SPEC, migclass, Foraging, propFor, volume)

write.csv(spp_table_traits, "traits/forest_spp_traits_MS.csv", row.names = F)

# Species traits plots
ed_breadth <- ggplot(env_breadth, aes(x = reorder(SPEC, ed), y = ed, color = SPEC)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ed - 1.96*std_ed, ymax = ed + 1.96*std_ed)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Forest edge density") +
  scale_color_viridis_d() +
  theme(legend.position = "none")

cover_breadth <- ggplot(env_breadth_allroutes, aes(x = reorder(SPEC, propFor), y = propFor)) +
  geom_point(cex = 8) + 
  scale_color_viridis_d()+
  geom_errorbar(aes(ymin = min_for, ymax = max_for), cex = 3, width = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Forest cover") +
  theme(legend.position = "none", axis.text.x = element_text(size = 26), 
        axis.text.y = element_text(size = 40),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 45),
        axis.line = element_line(colour = 'black', size = 3))

volume_plot <- ggplot(volume, aes(x = reorder(SPEC, volume), y = volume)) +
  geom_point(cex = 8) + 
  scale_color_viridis_d()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Climatic niche breadth") +
  theme(legend.position = "none", axis.text.x = element_text(size = 26), 
        axis.text.y = element_text(size = 40), 
        axis.title = element_text(size = 45),
        axis.line = element_line(colour = 'black', size = 3))

plot_grid(cover_breadth, volume_plot, nrow = 2)
ggsave("figures/main_analysis_figs/forest_breadth.tiff", units = "in", height = 13, width = 30)
ggsave("figures/main_analysis_figs/forest_breadth.pdf", units = "in", height = 13, width = 30)

## Correlations between forest edge density breadth, forest cover breadth, climate niche breadth

ed_cov <- ggplot(spp_breadths, aes(x = ed, y = propFor)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  labs(x = "Mean forest edge density", y = "Mean forest cover")

ed_vol <- ggplot(spp_breadths, aes(x = ed, y = volume)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Mean forest edge density", y = "Climatic niche breadth")

vol_cov <- ggplot(spp_breadths, aes(x = volume, y = propFor)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Climatic niche breadth", y = "Mean forest cover")

plot_grid(ed_cov, ed_vol, vol_cov, nrow = 2)
ggsave("figures/main_analysis_figs/traits_covariation.tiff", units = "in", height = 8, width = 8)

cor(spp_breadths[, c(4,7,9)]) # edxpropFor = -0.59, edxVolume = 0.08, volumexPropFor = -0.56

#### Individual models ####
# of climate + frag + loss + climate:frag + climate:loss for species

z <- function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T)}

forest_z <- data.frame(stateroute = forest$stateroute,
                       ED = z(forest$ED),
                       propForest = z(forest$propForest),
                       deltaED = z(forest$deltaED), 
                       deltaProp = z(forest$deltaProp))

climate_wide_z <- climate_wide %>%
  mutate_at(.vars = c("tmin", "tmax"),
            .funs = z)

clim_hab_poptrend_z <- abund_trend %>%
  left_join(dplyr::select(traits.short, -Common_name)) %>%
  left_join(forest_z, by = "stateroute") %>%
  left_join(climate_wide_z, by = "stateroute") %>%
  filter(!is.na(propForest)) %>%
  left_join(fourletter_codes)

# Spatial CAR models by species

model_fits <- clim_hab_poptrend_z %>%
  group_by(aou, SPEC) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .

    # get species weights matrix

    # Coordinates of BBS routes

    coords_bbs <- df %>%
      ungroup() %>%
      dplyr::select(stateroute) %>%
      distinct() %>%
      left_join(routes) %>%
      dplyr::select(stateroute, longitude, latitude)

    coords_mat <- as.matrix(coords_bbs[, -1])

    # Calculate nearest neighbors and make spatial neighborhood
    k0 <- knearneigh(coords_mat, longlat=T, k=1)
    k1 <- knn2nb(k0)

    # Find maximum neighbor distance and use this distance to define neighborhoods
    max.d <- max(unlist(nbdists(k1, coords_mat, longlat=T)))
    nb0.birds <- dnearneigh(coords_mat, 0, max.d, longlat=T)
    plot(nb0.birds, coords_mat)

    # Using a distance threshold of 100km
    nb1.birds <- dnearneigh(coords_mat,1,100,longlat=T)
    plot(nb1.birds, coords_mat)

    # Create spatial weights based on linear distance decay
    glist.birds <- nbdists(nb0.birds, coords_mat, longlat=T)
    glist.birds <- lapply(glist.birds, function(x) 1-(x/ceiling(max.d))) # Round up max.d so that all point get weights
    wt.birds <- nb2listw(nb0.birds, style='B', glist=glist.birds)

    spautolm(abundTrend ~ tmax*deltaED + tmin*deltaED + deltaProp, data = df, wt.birds, na.action = na.fail, family = "CAR")
  })) %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    nrow(df)
  }))  %>%
  mutate(lm_broom = map(lmFit,  ~{
    mod <- .
    sum <- summary(mod)
    df <- as.data.frame(sum$Coef)
    df$term = row.names(df)
    df
  })) %>%
  dplyr::select(aou, SPEC, nObs, lm_broom) %>%
  unnest() %>%
  filter(nObs > 40) %>% # 1 spp at 38, only one below 40
  filter(term != "(Intercept)") %>%
  rename(std.error = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`) %>%
  mutate(sig = case_when(p.value < 0.05 ~ "p < 0.05",
                         TRUE ~ "p > 0.05"))

range(model_fits$nObs)
# 67 species

# write.csv(model_fits, "model/individual_species_model_tables.csv", row.names = F)

model_fits <- read.csv("model/individual_species_model_tables.csv", stringsAsFactors = F)

# Figure: distributions of effect sizes by species

model_fits <- model_fits %>%
  mutate(sig = case_when(p.value < 0.05 & p.value >= 0.01 ~ "0.01 < p < 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "0.001 < p < 0.01",
                         p.value < 0.001 ~ "p < 0.001",
                         TRUE ~ "p > 0.05")) %>%
  mutate(sig = fct_relevel(sig, c("p < 0.001", "0.001 < p < 0.01", "0.01 < p < 0.05", "p > 0.05"))) %>%
  mutate(dir = ifelse(Estimate < 0, "Negative effect", "Positive effect"))

density_plot <- function(df, variable, label) {
  ggplot(filter(df, term == variable), aes(x = Estimate, fill = sig)) + 
    geom_histogram(bins = 20) + 
    geom_vline(aes(xintercept = mean(Estimate)), lty = 2, cex = 1) + 
    labs(x = label, y = "Species", fill = "") +
    scale_fill_manual(values = c("p < 0.001" = "#2166AC", 
                                 "0.001 < p < 0.01" = "#67A9CF",
                                 "0.01 < p < 0.05" = "#D1E5F0",
                                 "p > 0.05" = "gray"))
}

deltaED <- density_plot(model_fits, "deltaED", "Change in edge density")
deltaProp <- density_plot(model_fits, "deltaProp", "Change in forest cover") + scale_x_continuous(breaks = c(-0.01, 0, 0.01))
tmin <- density_plot(model_fits, "tmin", "Trend in Tmin")
tmax <- density_plot(model_fits, "tmax", "Trend in Tmax")
tmaxED <- density_plot(model_fits, "tmax:deltaED", "Tmax:change in ED")
tminED <- density_plot(model_fits, "deltaED:tmin", "Tmin:change in ED")

legend <- get_legend(deltaED)

grid_effects <- plot_grid(deltaED + theme(legend.position = "none"), 
                          deltaProp + theme(legend.position = "none") + ylab(" "), 
          tmin + theme(legend.position = "none"), 
          tmax + theme(legend.position = "none") + ylab(" "), 
          tmaxED + theme(legend.position = "none") + scale_x_continuous(breaks = c(-0.01, 0, 0.01)), 
          tminED + theme(legend.position = "none") + ylab(" "),
          nrow = 3,
          labels = c("A", "B", "C", "D", "E", "F"))
plot_grid(grid_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/model_effects_distributions.pdf", units = "in", height = 9, width = 10)

## Supplemental figure: p value ranges for each variable

pval_plot <- function(df, variable, label) {
  ggplot(filter(df, term == variable), aes(x = p.value, fill = dir)) + 
    geom_histogram(bins = 20) + 
    geom_vline(aes(xintercept = mean(p.value)), lty = 2, cex = 1) + 
    labs(x = label, y = "Count", fill = "")  +
    scale_x_log10() + scale_fill_manual(values = c("dodgerblue3", "firebrick2"))

}

options(scipen = 999)
deltaED <- pval_plot(model_fits, "deltaED", "Change in edge density") + scale_x_log10(breaks = c(0.00001, 0.001, 1),
                                                                labels = c(0.00001, 0.001, 1))
deltaProp <- pval_plot(model_fits, "deltaProp", "Change in forest cover")
tmin <- pval_plot(model_fits, "tmin", "Trend in Tmin")
tmax <- pval_plot(model_fits, "tmax", "Trend in Tmax")
tmaxED <- pval_plot(model_fits, "tmax:deltaED", "Tmax:change in ED")
tminED <- pval_plot(model_fits, "deltaED:tmin", "Tmin:change in ED")

legend <- get_legend(tmin)

grid_effects <- plot_grid(deltaED + theme(legend.position = "none"), 
                          deltaProp + ylab(" ") + theme(legend.position = "none"), 
                          tmin + theme(legend.position = "none"), 
                          tmax + ylab(" ") + theme(legend.position = "none"), 
                          tmaxED + theme(legend.position = "none"), 
                          tminED + ylab(" ") + theme(legend.position = "none"),
                          nrow = 3,
                          labels = c("A", "B", "C", "D", "E", "F"))
plot_grid(grid_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/model_pvals_distributions.pdf", units = "in", height = 9, width = 10)

## Traits and responses 

forest_traits <- spp_breadths %>%
  left_join(area_aous) %>%
  left_join(correlates, by = c("aou" = "AOU"))

spp_traits <- spp_breadths %>%
  left_join(area_aous) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  left_join(model_fits) %>%
  group_by(term) %>%
  nest() %>%
  filter(!is.na(term)) %>%
  mutate(trait_mod = purrr::map(data, ~{
    df <- .
    lm(Estimate ~  volume + propFor + migclass + Foraging, data = df)
  }),
  tidy = purrr::map(trait_mod, ~{
    mod <- .
    tidy(mod)
  }),
  r2 = map_dbl(trait_mod, ~{
    mod <- .
    glance(mod)$r.squared
  })) %>%
  dplyr::select(term, tidy, r2) %>%
  unnest()

# write.csv(spp_traits, "model/spp_trait_model_output.csv", row.names = F)
spp_traits <- read.csv("model/spp_trait_model_output.csv", stringsAsFactors = F)

spp_traits_pred <- spp_traits %>%
  filter(p.value < 0.05)

### species traits models effect plots

trait_effects <- spp_breadths %>%
  left_join(area_aous) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  left_join(model_fits)

# propFor v deltaED
ggplot(filter(trait_effects, term == "deltaED"), aes(x = propFor, y = Estimate)) + 
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(y = "Effect of change in forest edge density", x = "Mean proportion forest cover")
ggsave("figures/main_analysis_figs/trait_model_deltaED_propFor.pdf")

# volume v tmin

ggplot(filter(trait_effects, term == "tmin"), aes(x = volume, y = Estimate)) + 
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(y = "Effect of trend in Tmin", x = "Climatic niche breadth")
ggsave("figures/main_analysis_figs/trait_model_tmin_volume.pdf")

# Range centroid for each species - map of effect sizes for each predictor

species_centroids <- counts.subs %>%
  ungroup() %>%
  filter(aou %in% model_fits$aou) %>%
  group_by(aou) %>%
  distinct(stateroute, latitude, longitude) %>%
  summarize(centroid_lon = mean(longitude), centroid_lat = mean(latitude)) %>%
  st_as_sf(coords = c("centroid_lon", "centroid_lat"))

terms <- unique(model_fits$term)[1:4]

for(trm in terms) {
  df <- model_fits %>%
    filter(term == trm)
  
  species_centroids_efx <- species_centroids %>%
    left_join(df, by = "aou") %>%
    mutate(dot_size = abs(Estimate))
  
  map <- tm_shape(na_map) + tm_fill() + tm_shape(species_centroids_efx) + 
    tm_dots(col = "dir", size = "dot_size", palette = "RdBu", scale = 2, alpha = 0.7) +
    tm_layout(main.title = trm)
  tmap_save(map, paste0("figures/", trm, "_effects_map_spp_centroids.pdf"))
}

#### Range position models ####

## Mean breeding season temperature on BBS routes

route_climate <- read.csv("climate/bbs_routes_breeding_season_climate.csv", stringsAsFactors = F) %>%
  group_by(stateroute) %>%
  summarize(mean_tmax = mean(mean_tmax),
            mean_tmin = mean(mean_tmin))

## Join route position with clim_hab_poptrend table

clim_hab_poptrend_means <- clim_hab_poptrend_z %>%
  left_join(route_climate)

## Individual models of Tmin and Tmax with interactions with breeding season mean Tmin and Tmax over study period

 model_fits_position <- clim_hab_poptrend_means %>%
   group_by(aou, SPEC) %>%
   nest() %>%
   mutate(lmFit = map(data, ~{
     df <- .

     # get species weights matrix

     # Coordinates of BBS routes

     coords_bbs <- df %>%
       ungroup() %>%
       dplyr::select(stateroute) %>%
       distinct() %>%
       left_join(routes) %>%
       dplyr::select(stateroute, longitude, latitude)

     coords_mat <- as.matrix(coords_bbs[, -1])

     # Calculate nearest neighbors and make spatial neighborhood
     k0 <- knearneigh(coords_mat, longlat=T, k=1)
     k1 <- knn2nb(k0)

     # Find maximum neighbor distance and use this distance to define neighborhoods
     max.d <- max(unlist(nbdists(k1, coords_mat, longlat=T)))
     nb0.birds <- dnearneigh(coords_mat, 0, max.d, longlat=T)
     plot(nb0.birds, coords_mat)

     # Using a distance threshold of 100km
     nb1.birds <- dnearneigh(coords_mat,1,100,longlat=T)
     plot(nb1.birds, coords_mat)

     # Create spatial weights based on linear distance decay
     glist.birds <- nbdists(nb0.birds, coords_mat, longlat=T)
     glist.birds <- lapply(glist.birds, function(x) 1-(x/ceiling(max.d))) # Round up max.d so that all point get weights
     wt.birds <- nb2listw(nb0.birds, style='B', glist=glist.birds)

     spautolm(abundTrend ~ tmax*mean_tmax + tmin*mean_tmin, data = df, wt.birds, na.action = na.fail, family = "CAR")
   })) %>%
   mutate(nObs = map_dbl(data, ~{
     df <- .
     nrow(df)
   }))  %>%
   mutate(lm_broom = map(lmFit,  ~{
     mod <- .
     sum <- summary(mod)
     df <- as.data.frame(sum$Coef)
     df$term = row.names(df)
     df
   })) %>%
   dplyr::select(aou, SPEC, nObs, lm_broom) %>%
   unnest() %>%
   rename(std.error = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`)

# write.csv(model_fits_position, "model/range_position_model_tables.csv", row.names = F)
model_fits_position <- read.csv("model/range_position_model_tables.csv", stringsAsFactors = F)

model_fits_position_fig <- model_fits_position %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(p.value < 0.05 & p.value >= 0.01 ~ "0.01 < p < 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "0.001 < p < 0.01",
                         p.value < 0.001 ~ "p < 0.001",
                         TRUE ~ "p > 0.05")) %>%
  mutate(sig = fct_relevel(sig, c("p < 0.001", "0.001 < p < 0.01", "0.01 < p < 0.05", "p > 0.05"))) %>%
  mutate(dir = ifelse(Estimate > 0, "Positive effect", "Negative effect")) %>%
  ungroup() %>%
  filter(SPEC != "CERW")

# Figure

tmax_int_plot <- density_plot(model_fits_position_fig, "tmax:mean_tmax", "Interaction: trend in Tmax & mean Tmax")
tmin_int_plot <- density_plot(model_fits_position_fig, "tmin:mean_tmin", "Interaction: trend in Tmin & mean Tmin")
  
legend <- get_legend(tmax_int_plot)

temp_effects <- plot_grid(tmax_int_plot + theme(legend.position = "none"), 
                          tmin_int_plot + theme(legend.position = "none") + ylab(" "),
                          nrow = 1,
                          labels = c("A", "B"))
plot_grid(temp_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/range_position_temp_responses.pdf", units = "in", height = 3, width = 10)

## Supplemental figure: p-values

tmax_int_plot <- pval_plot(model_fits_position_fig, "tmax:mean_tmax", "Interaction: trend in Tmax & mean Tmax")
tmin_int_plot <- pval_plot(model_fits_position_fig, "tmin:mean_tmin", "Interaction: trend in Tmin & mean Tmin")

legend <- get_legend(tmax_int_plot)

temp_effects <- plot_grid(tmax_int_plot + theme(legend.position = "none"), 
                          tmin_int_plot + theme(legend.position = "none") + ylab(" "),
                          nrow = 1,
                          labels = c("A", "B"))
plot_grid(temp_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/range_position_temp_pvals.pdf", units = "in", height = 3, width = 10)

#### Strongest species response for each predictor ####

# NOFL - deltaED

flicker <-  clim_hab_poptrend_z %>%
  filter(Common_name == "Northern flicker")

flicker_plot <- ggplot(flicker, aes(x = deltaED, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Change in edge density", y = "Abundance trend", title = "Northern flicker")

# INBU - deltaProp

bunting <- clim_hab_poptrend_z %>%
  filter(Common_name == "Indigo bunting") 

bunting_plot <- ggplot(bunting, aes(x = deltaProp, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Change in forest cover", y = "", title = "Indigo bunting")
 
# BAWW - tmax

warbler <- clim_hab_poptrend_z %>%
  filter(Common_name == "Black-and-white warbler")

warbler_plot <- ggplot(warbler, aes(x = tmax, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Trend in Tmax", y = "", title = "Black-and-white warbler")

# MAWA - tmin

magnolia <- clim_hab_poptrend_z %>%
  filter(Common_name == "Magnolia warbler")

magnolia_plot <- ggplot(magnolia, aes(x = tmin, y = abundTrend)) + geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Trend in Tmin", y = "Abundance trend", title = "Magnolia warbler")

# SUTA - deltaED:tmax

tanager <- clim_hab_poptrend_z %>%
  filter(Common_name == "Summer tanager") 

tanager$tmax_sign <- ifelse(tanager$tmax < 0, "-1.5 < Z-Tmax < 0", "0 < Z-Tmax < 1.7")
tanager_plot <- ggplot(tanager, aes(x = deltaED, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~tmax_sign) +
  theme(panel.spacing = unit(4, "lines")) +
  labs(x = "Change in ED", y = "Abundance trend", title = "Summer tanager")

## Figure for MS

indiv_spp_add <- plot_grid(flicker_plot, bunting_plot, warbler_plot, magnolia_plot, nrow = 2, labels = c("A", "B", "C", "D"))
indiv_spp_multi <- plot_grid(indiv_spp_add, tanager_plot, nrow = 2, rel_heights = c(0.66, 0.33),
                             labels = c(" ", "E", "F"))
ggsave("figures/main_analysis_figs/indiv_spp_multipanel.pdf", indiv_spp_multi, units = "in", width = 9, height = 10)

