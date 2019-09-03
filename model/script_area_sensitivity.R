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
library(MuMIn)
library(nlme)
library(grid)
library(spdep)

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
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# species four letter codes
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
spp_codes <- read.csv("four_letter_codes_birdspp.csv", stringsAsFactors = F)
spp_codes$common_name_lower <- tolower(spp_codes$COMMONNAME)

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)
setwd("\\\\BioArk/hurlbertlab/Databases/BBS/GPS_stoplocations/")
setwd("/Volumes/hurlbertlab/Databases/BBS/GPS_stoplocations/")
us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data$rteno > 3000), ]

# Subset stateroutes: runtype = 1, sampled once every 5 years since 1990, between 38 and 42 km in length
routes.short <- RT1.routes %>%
  filter(stateroute %in% us_subs@data$rteno) %>%
  mutate(year_bin = 5*floor(year/5) + 5/2) %>%
  filter(year >= 1990) %>%
  group_by(stateroute) %>%
  distinct(year_bin) %>%
  count() %>%
  filter(n >= 5)
# 1513 routes

# Subset species: forest breeding birds by area sensitivity

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
  filter(year >= 1990, year < 2017)
# 1513 routes

# first year observers

first_year <- function(obsns) {
  first_yrs <- c()
  for(i in 1:length(obsns)) {
    ob <- obsns[i]
    obs <- obsns[1:i]
    first_yrs <- c(first_yrs, ifelse(length(obs[obs == ob]) == 1, 1, 0))
  }
  return(first_yrs)
}

obs_years <- weather %>%
  group_by(stateroute) %>%
  arrange(year) %>%
  nest() %>%
  mutate(first_yr = purrr::map(data, ~{
    df <- .
    first_year(df$obsn)
  })) %>%
  unnest() %>%
  dplyr::select(stateroute, year, first_yr)

# abundance trends read in

abund_trend <- counts.subs %>%
  left_join(obs_years, by = c('stateroute', 'year')) %>%
  group_by(aou, stateroute) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    df.short <- df %>%
      dplyr::select(year, first_yr, speciestotal) %>%
      unique()
    glm(speciestotal ~ year + first_yr, family = poisson, data = df.short)
  })) %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    nrow(df)
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  mutate(abundTrend = map_dbl(lm_broom, ~{
    df <- .
    df$estimate[2]
  })) %>%
  mutate(trendInt = map_dbl(lm_broom, ~{
    df <- .
    df$estimate[1]
  })) %>%
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    df$p.value[2]
  })) %>%
  mutate(obsTrend = map_dbl(lm_broom, ~{
    df <- .
    df$estimate[3]
  })) %>%
  mutate(obsPval = map_dbl(lm_broom, ~{
    df <- .
    df$p.value[3]
  }))

hist(abund_trend$abundTrend)
hist(abund_trend$nObs)

# Save abundance trend data
abund_trend <- abund_trend %>%
  filter(nObs > 9) %>%
  dplyr::select(-data, -lmFit, -lm_broom)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
write.csv(abund_trend, "BBS_abundance_trends.csv", row.names = F)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
abund_trend <- read.csv("BBS_abundance_trends.csv", stringsAsFactors = F) %>%
  filter(aou %in% area_aous$aou) %>% ## Only ones lost from area_aous are nocturnal: owls, nightjars (Caprimulgus)
  left_join(area_aous) %>%
  filter(sporder != "Accipitriformes" & sporder != "Strigiformes" & sporder != "Caprimulgiformes") ## Remove hawks

length(abund_trend$trendPval[abund_trend$trendPval < 0.05])/nrow(abund_trend) # 44%
length(abund_trend$trendPval[abund_trend$trendPval < 0.01])/nrow(abund_trend) # 35%
length(abund_trend$trendPval[abund_trend$trendPval < 0.0001])/nrow(abund_trend) # 21%

abund_trend_sig <- abund_trend %>%
  mutate(abund_class = case_when(trendPval >= 0.05 ~ "p > 0.05",
                                 trendPval >= 0.01 & trendPval < 0.05 ~ "0.01 < p < 0.05",
                                 trendPval >= 0.001 & trendPval < 0.01 ~ "0.001 < p < 0.01",
                                 trendPval < 0.001 ~ "p < 0.001"),
         obs_class = case_when(obsPval >= 0.05 ~ "p > 0.05",
                               obsPval >= 0.01 & obsPval < 0.05 ~ "0.01 < p < 0.05",
                               obsPval >= 0.001 & obsPval < 0.01 ~ "0.001 < p < 0.01",
                               obsPval < 0.001 ~ "p < 0.001"))

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")

# Histogram of abundance trends colored by p value significance level

ggplot(abund_trend_sig, aes(x = abundTrend, 
                            fill = fct_relevel(abund_class, c("p < 0.001",
                                                              "0.001 < p < 0.01", 
                                                              "0.01 < p < 0.05", 
                                                              "p > 0.05")))) + 
  geom_histogram(bins = 20) + 
  labs(fill = "", x = "Abundance trend", y = "Count") + 
  scale_fill_viridis_d()
ggsave("figures/area_sensitivity/hist_abundTrend_sig.pdf")

# Individual maps for species of abundance trends in range

us.proj <- readOGR("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/BCRs_contiguous_us.shp")
states <- us.proj[, -(1:2)]

abund_trend_map <- abund_trend %>%
  left_join(routes) %>%
  st_as_sf(coords = c("longitude", "latitude"))

pdf("figures/area_sensitivity/species_abundance_maps.pdf")
for(spp in unique(abund_trend$aou)) {
  df <- filter(abund_trend_map, aou == spp)
  spec <- fourletter_codes$SPEC[fourletter_codes$aou == spp]

  print(tm_shape(states) + tm_polygons(col = 'white') + tm_shape(df) + tm_dots(col = "abundTrend", palette = "RdBu", size = 1) +
    tm_layout(title = spec))
}
dev.off()

# Individual species maps - average by state

states_sf <- read_sf("\\\\BioArk\\Hurlbertlab\\GIS\\geography\\ne_50m_admin_1_states_provinces_lakes\\ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA") %>%
  dplyr::select(adm1_code) %>%
  filter(adm1_code != "USA-3563" & adm1_code != "USA-3517")
crs <- st_crs(states_sf)

pdf("figures/area_sensitivity/species_abundance_maps_states.pdf")
for(spp in unique(abund_trend$aou)) {

  spec <- fourletter_codes$SPEC[fourletter_codes$aou == spp]
  
  abund_states <- abund_trend_map %>%
    filter(aou == spp) %>%
    st_set_crs(crs) %>%
    st_join(states_sf) %>%
    group_by(adm1_code) %>%
    summarize(meanTrend = mean(abundTrend, na.rm = T))
  
  abund_states_df <- data.frame(adm1_code = abund_states$adm1_code, meanTrend = abund_states$meanTrend)
  
  abund_states_polygon <- states_sf %>%
    left_join(abund_states_df)
  
  print(tm_shape(abund_states_polygon) + tm_borders() + 
          tm_fill(col= "meanTrend", palette = "RdBu") + tm_layout(title = spec))
 
}
dev.off()

# For one species with not that many routes - show counts vs. time

bhvi <- filter(abund_trend, aou == 6290, nObs > 9)

pdf("figures/area_sensitivity/abundance_trends_bhvi.pdf")
for(strte in bhvi$stateroute) {
  df <- filter(bhvi, stateroute == strte)
  glm <- df$lmFit[[1]]
  df.count <- df$data[[1]]
  df.count$predict <- predict(glm, newdata = df.count, type = "response")
  plot <- ggplot(df.count, aes(x = year, y = speciestotal)) + geom_point() + 
    geom_smooth(method = "glm", se = F, method.args = list(family = "poisson"))
    
  print(plot)

  }
dev.off()

# Histogram of first year observer effect colored by p value significance level

ggplot(abund_trend_sig, aes(x = obsTrend, 
                            fill = fct_relevel(obs_class, c("p < 0.001",
                                                              "0.001 < p < 0.01", 
                                                              "0.01 < p < 0.05", 
                                                              "p > 0.05")))) + 
  geom_histogram(bins = 20) + 
  labs(fill = "", x = "First year observer effect", y = "Count") + 
  scale_fill_viridis_d()
ggsave("figures/area_sensitivity/hist_obsEffect_sig.pdf")

# Which species most commonly show observer effects?

presences <- counts.subs %>%
  group_by(aou) %>%
  summarize(range_size = n_distinct(stateroute))

obs_effects <- abund_trend_sig %>% 
  filter(!is.na(obs_class), obsPval < 0.05) %>%
  group_by(aou) %>%
  count() %>%
  left_join(presences) %>%
  left_join(fourletter_codes) %>%
  mutate(obs_prop = n/range_size) %>%
  arrange(desc(obs_prop))
write.csv(obs_effects, "model/obs_effect_rankings.csv", row.names = F)

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

# Landcover legend
newcode <- data.frame(code = seq(1,9), 
                      legend = c("Open water", "Urban", "Barren", "Forest", "Shrubland", 
                                 "Agricultural", "Grasslands", "Wetlands", "Perennial ice, snow"))

# Forest fragmentation deltas
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

forest <- forest_ed %>%
  right_join(forest_deltaED) %>%
  left_join(forest_deltaP)

# Master data table
clim_hab_poptrend <- abund_trend %>%
  left_join(dplyr::select(traits.short, -Common_name)) %>%
  left_join(forest, by = "stateroute") %>%
  left_join(climate_wide, by = "stateroute") %>%
  filter(!is.na(propForest)) %>%
  group_by(aou) %>% 
  mutate(abundTrend_z = (abundTrend - mean(abundTrend, na.rm = T))/sd(abundTrend, na.rm = T))

##### Analysis #####
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/")

#### Route trends ####

# Plot of routes used

forest_routes <- us_subs[us_subs$rteno %in% forest_ed$stateroute, ]

us.proj <- readOGR("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/BCRs_contiguous_us.shp")

us_subs_transf <- spTransform(forest_routes, crs(us.proj))

plot(us.proj, col = "gray73", border = "gray73")
plot(us_subs_transf, add = T)


hist(clim_hab_poptrend$deltaED)
hist(clim_hab_poptrend$deltaProp)

cor(clim_hab_poptrend$ED, clim_hab_poptrend$propForest) # -0.31
cor(clim_hab_poptrend$deltaED, clim_hab_poptrend$deltaProp) # -0.24
cor(dplyr::select(clim_hab_poptrend, ppt, tmin, tmax, deltaED, deltaProp, ED, propForest))

# Figure: map of routes

color_scale <- data.frame(color = c(1:5), 
                          dED_hex = c("#276419", "#7fbc41", "#e6f5d0", "#fde0ef", "#de77ae"),
                          dProp_hex = c("#7b3294", "#c2a5cf", "#e7d4e8", "#d9f0d3", "#5aae61"),
                          temp_hex = c("#2166AC", "#7aafcf", "#fddbc7", "#f4a582", "#b2182b"), stringsAsFactors = F)

states <- us.proj[, -(1:2)]

bbs_routes_forest <- us_subs_transf %>%
  st_as_sf() %>%
  left_join(forest, by = c("rteno" = "stateroute")) %>%
  left_join(climate_wide, by = c("rteno" = "stateroute")) %>%
  filter(!is.na(deltaProp))

bbs_routes_forcov <- tm_shape(states) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + 
  tm_lines(col = "deltaProp", scale = 3, palette = "PRGn") +
  tm_layout(title = "A", title.fontface = "bold", legend.show = F)

bbs_routes_fored <- tm_shape(states) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + tm_lines(col = "deltaED", scale = 3, palette = "-PiYG") +
  tm_layout(legend.show = F, title = "B", title.fontface = "bold")

bbs_routes_tmin <- tm_shape(states) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + tm_lines(col = "tmin", scale = 3, breaks = c(-0.07, -0.05, 0, 0.05, 0.1, 0.2), palette = color_scale$temp_hex) +
  tm_layout(legend.show = F, title = "C", title.fontface = "bold")

bbs_routes_tmax <- tm_shape(states) + tm_fill(col = "gray63") + tm_borders(col = "gray63") + 
  tm_shape(bbs_routes_forest)  + tm_lines(col = "tmax", scale = 3, palette = color_scale$temp_hex) +
  tm_layout(legend.show = F, title = "D", title.fontface = "bold")

bbs_routes <- tmap_arrange(bbs_routes_forcov, bbs_routes_fored, bbs_routes_tmin, bbs_routes_tmax, nrow = 2)

# Histogram legends for each variable

# Add color scale

bbs_routes_forest_plots <- bbs_routes_forest %>%
  mutate(dED_color = case_when(
    deltaED <= -0.2 ~ 1,
    deltaED > -0.2 & deltaED <= -0.1 ~ 2,
    deltaED > -0.1 & deltaED <= 0 ~ 3,
    deltaED > 0 & deltaED <= 0.1 ~ 4,
    deltaED > 0.1 ~ 5
  ),
  dProp_color = case_when(
    deltaProp <= -0.4 ~ 1,
    deltaProp > -0.4 & deltaProp <= -0.2 ~ 2,
    deltaProp > -0.2 & deltaProp <= 0 ~ 3,
    deltaProp > 0 & deltaProp <= 0.2 ~ 4,
    deltaProp > 0.2 ~ 5
  ),
  tmin_color = case_when(
    tmin <= -0.05 ~ 1,
    tmin > -0.05 & tmin <= 0 ~ 2,
    tmin > 0 & tmin <= 0.05 ~ 3,
    tmin > 0.05 & tmin <= 0.1 ~ 4,
    tmin > 0.1 ~ 5
  ),
  tmax_color = case_when(
    tmax <= -0.05 ~ 1,
    tmax > -0.05 & tmax <= 0 ~ 2,
    tmax > 0 & tmax <= 0.05 ~ 3,
    tmax > 0.05 & tmax <= 0.1 ~ 4,
    tmax > 0.1 ~ 5
  ))

ED_hist <- ggplot(bbs_routes_forest_plots, aes(x = deltaED, fill = as.factor(dED_color))) + geom_histogram(bins = 50) +
  labs(title = "Change in edge density") +
  scale_y_continuous(breaks = c(0, 75, 150)) +
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
  dplyr::select(-meanPatchArea, -stateroute, -ppt, -ED, -propForest, -year, -deltaEDrecent, -recentChange) %>%
  rename("Tmin" = "tmin", "Tmax" = "tmax", "forestED" = "deltaED", "forestCover" = "deltaProp")

corr.table <- cor(matrix)
pdf("figures/route_eda/env_cor_matrix.pdf")
corrplot::corrplot(corr.table, method = "circle", diag = F, tl.col = "black",
                   tl.cex = 1.5, cl.cex = 1)
dev.off()

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

forest_deltaP_allroutes <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(propForest = prop.landscape) %>%
  spread(key = "year", value = "propForest") %>%
  group_by(stateroute) %>%
  summarize(deltaProp = `2016` - `1992`)

forest_allroutes <- forest_ed_allroutes %>%
  left_join(forest_deltaP_allroutes)

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
ggsave("figures/area_sensitivity/forest_breadth.tiff", units = "in", height = 13, width = 30)
ggsave("figures/area_sensitivity/forest_breadth.pdf", units = "in", height = 13, width = 30)

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
ggsave("figures/area_sensitivity/traits_covariation.tiff", units = "in", height = 8, width = 8)

cor(spp_breadths[, c(4,7,9)]) # edxpropFor = -0.61, edxVolume = 0.13, volumexPropFor = -0.60

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
  group_by(SPEC, aou) %>% 
  mutate(abundTrend_z = (abundTrend - mean(abundTrend, na.rm = T))/sd(abundTrend, na.rm = T)) %>%
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
  filter(nObs > 40) %>%
  filter(term != "(Intercept)") %>%
  rename(std.error = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`) %>%
  mutate(sig = case_when(p.value < 0.05 ~ "p < 0.05",
                         TRUE ~ "p > 0.05"))

range(model_fits$nObs) # 1 spp at 38, only one below 40
# 67 species

#write.csv(model_fits, "model/individual_species_model_tables.csv", row.names = F)

model_fits <- read.csv("model/individual_species_model_tables.csv", stringsAsFactors = F)

# Compare CAR model with linear model results

model_fits_linear <- clim_hab_poptrend_z %>%	
  group_by(aou, SPEC) %>%	
  nest() %>%	 
  mutate(lmFit = purrr::map(data, ~{
    df <- .	   
    lm(abundTrend ~ tmax*deltaED + tmin*deltaED + deltaProp, df, na.action = na.fail)
  })) %>%	
  mutate(nObs = map_dbl(data, ~{	
    df <- .	  
    nrow(df)
  })) %>%
  mutate(lm_broom = purrr::map(lmFit,  ~{
    tidy(.)
  })) %>%
  dplyr::select(aou, SPEC, nObs, lm_broom) %>%	 
  unnest() %>%
  filter(nObs > 40) %>%	 
  filter(term != "(Intercept)")

# Join with CAR

for(trm in unique(model_fits_linear$term)) {
  name <- str_remove(trm, ":")
  linear <- filter(model_fits_linear, term == trm) %>%
    mutate(model = "linear") %>%
    dplyr::select(aou, SPEC, estimate, model)
  car <- filter(model_fits, term == trm) %>%
    mutate(model = "car") %>%
    dplyr::select(aou, SPEC, Estimate, model) %>%
    rename(estimate = Estimate)
  
  mods <- bind_rows(linear, car) %>%
    spread(model, estimate)
  
  plot <- ggplot(mods, aes(x = car, y = linear)) + geom_point() + geom_abline(slope = 1) + 
    geom_text(aes(label = SPEC), nudge_y = 0.0005)
  ggsave(paste0("figures/area_sensitivity/", name, "_comparison.pdf"), plot)
}

# Figure: distributions of effect sizes by species

model_fits <- model_fits %>%
  mutate(sig = case_when(p.value < 0.05 & p.value >= 0.01 ~ "0.01 < p < 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "0.001 < p < 0.01",
                         p.value < 0.001 ~ "p < 0.001",
                         TRUE ~ "p > 0.05")) %>%
  mutate(sig = fct_relevel(sig, c("p < 0.001", "0.001 < p < 0.01", "0.01 < p < 0.05", "p > 0.05"))) %>%
  mutate(dir = ifelse(Estimate < 0, "Negative effect", "Positive effect"))

density_plot <- function(variable, label) {
  ggplot(filter(model_fits, term == variable), aes(x = Estimate, fill = sig)) + 
    geom_histogram(bins = 20) + 
    geom_vline(aes(xintercept = mean(Estimate)), lty = 2, cex = 1) + 
    labs(x = label, y = "Species", fill = "") +
    scale_fill_manual(values = c("p < 0.001" = "#2166AC", 
                                 "0.001 < p < 0.01" = "#67A9CF",
                                 "0.01 < p < 0.05" = "#D1E5F0",
                                 "p > 0.05" = "gray"))
}

deltaED <- density_plot("deltaED", "Change in ED")
deltaProp <- density_plot("deltaProp", "Change in forest cover") + scale_y_continuous(breaks = c(0,5,10))
tmin <- density_plot("tmin", "Trend in Tmin")
tmax <- density_plot("tmax", "Trend in Tmax") + scale_y_continuous(breaks = c(0, 4, 8))
tmaxED <- density_plot("tmax:deltaED", "Tmax:change in ED")
tminED <- density_plot("deltaED:tmin", "Tmin:change in ED")

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
ggsave("figures/area_sensitivity/model_effects_distributions.pdf", units = "in", height = 9, width = 10)

## Supplemental figure: p value ranges for each variable

pval_plot <- function(variable, label) {
  ggplot(filter(model_fits, term == variable), aes(x = p.value, fill = dir)) + 
    geom_histogram(bins = 20) + 
    geom_vline(aes(xintercept = mean(p.value)), lty = 2, cex = 1) + 
    labs(x = label, y = "Count", fill = "")  +
    scale_x_log10() + scale_fill_manual(values = c("dodgerblue3", "firebrick2"))

}

options(scipen = 999)
deltaED <- pval_plot("deltaED", "Change in ED") + scale_x_log10(breaks = c(0.00001, 0.001, 1),
                                                                labels = c(0.00001, 0.001, 1))
deltaProp <- pval_plot("deltaProp", "Change in forest cover") + scale_x_log10(breaks = c(0.0001, 0.01, 1),
                                                                             labels = c(0.0001, 0.01, 1))
tmin <- pval_plot("tmin", "Trend in Tmin") + scale_x_log10(breaks = c(0.0001, 0.01, 1),
                                                           labels = c(0.0001, 0.01, 1))
tmax <- pval_plot("tmax", "Trend in Tmax")
tmaxED <- pval_plot("tmax:deltaED", "Tmax:change in ED")
tminED <- pval_plot("deltaED:tmin", "Tmin:change in ED")

grid_effects <- plot_grid(deltaED, 
                          deltaProp + ylab(" "), 
                          tmin, 
                          tmax + ylab(" "), 
                          tmaxED, 
                          tminED + ylab(" "),
                          nrow = 3,
                          labels = c("A", "B", "C", "D", "E", "F"))
plot_grid(grid_effects)
ggsave("figures/area_sensitivity/model_pvals_distributions.pdf", units = "in", height = 9, width = 10)

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

#write.csv(spp_traits, "model/spp_trait_model_output.csv", row.names = F)
spp_traits <- read.csv("model/spp_trait_model_output.csv", stringsAsFactors = F)

spp_traits_pred <- spp_traits %>%
  filter(p.value < 0.05)

### species traits models effect plots

trait_effects <- spp_breadths %>%
  left_join(area_aous) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  left_join(model_fits)

# volume v tmax
ggplot(filter(trait_effects, term == "tmax"), aes(x = volume, y = Estimate)) + 
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(y = "Effect of trend in Tmax", x = "Niche volume")
ggsave("figures/area_sensitivity/trait_model_tmax_volume.pdf")

# propFor v tmin
ggplot(filter(trait_effects, term == "tmin"), aes(x = propFor, y = Estimate)) + 
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(y = "Effect of trend in Tmin", x = "Mean breeding range forest cover")
ggsave("figures/area_sensitivity/trait_model_tmin_propFor.pdf")

#### Range position models ####

range_centroids <- read.csv("traits/species_Brange_centroids.csv", stringsAsFactors = F)
route_directions <- read.csv('traits/species_route_dist_from_centroids.csv', stringsAsFactors = F)

# Group routes into north, south, center of range

route_terciles <- route_directions %>%
  left_join(range_centroids, by = c("aou")) %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  mutate(range_direction = case_when(latitude >= n_tercile ~ "North", 
                                    latitude <= s_tercile ~ "South",
                                     TRUE ~ "Central")) %>%
  group_by(aou, range_direction) %>%
  filter(n() > 40)

## Join route position with clim_hab_poptrend table

clim_hab_poptrend_position <- clim_hab_poptrend_z %>%
  filter(aou %in% route_terciles$aou) %>%
  left_join(route_terciles, by = c("stateroute", "aou")) %>%
  filter(!is.na(range_direction)) %>%
  group_by(aou, range_direction) %>%
  mutate(n = n()) %>%
  filter(n > 40)

## Individual models of Tmin and Tmax with interactions with route position

model_fits_position <- clim_hab_poptrend_position %>%
  group_by(aou, SPEC, range_direction) %>%
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
    
    spautolm(abundTrend ~ tmax + tmin, data = df, wt.birds, na.action = na.fail, family = "CAR")
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
  dplyr::select(aou, SPEC, range_direction, nObs, lm_broom) %>%
  unnest() %>%
  rename(std.error = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`) 

write.csv(model_fits_position, "model/range_position_model_tables.csv", row.names = F)
model_fits_position <- read.csv("model/range_position_model_tables.csv", stringsAsFactors = F)

model_fits_position_fig <- model_fits_position %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(p.value < 0.05 & p.value >= 0.01 ~ "0.01 < p < 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "0.001 < p < 0.01",
                         p.value < 0.001 ~ "p < 0.001",
                         TRUE ~ "p > 0.05")) %>%
  mutate(sig = fct_relevel(sig, c("p < 0.001", "0.001 < p < 0.01", "0.01 < p < 0.05", "p > 0.05"))) %>%
  mutate(dir = ifelse(Estimate > 0, "Positive effect", "Negative effect")) %>%
  mutate(term = fct_recode(term, Tmin = "tmin", Tmax = "tmax")) %>%
  mutate(term = paste("Trend in", term)) %>%
  mutate(range_direction = fct_relevel(range_direction, levels = c("North", "Central", "South"))) %>%
  filter(range_direction != "Central") %>%
  ungroup() %>%
  filter(SPEC != "CERW") %>%
  group_by(SPEC) %>%
  filter(n_distinct(range_direction) == 2) %>%
  group_by(range_direction, term) %>%
  mutate(meanEst = mean(Estimate, na.rm =T))

# Figure

ggplot(model_fits_position_fig, aes(Estimate, fill = sig)) +
  geom_histogram(bins = 15) + 
  facet_grid(range_direction ~ term) +
  geom_vline(aes(xintercept = meanEst), lty = 2, cex = 1) +
  scale_fill_manual(values = c("p < 0.001" = "#2166AC", 
                               "0.001 < p < 0.01" = "#67A9CF",
                               "0.01 < p < 0.05" = "#D1E5F0",
                               "p > 0.05" = "gray")) +
  scale_x_continuous(breaks = c(0, 5, 10)) +
  labs(x = "Effect estimate",
       y = "Species",
       fill = "") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black")) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        text = element_text(size = 12),
        legend.text = element_text(size = 12))
  
ggsave("figures/area_sensitivity/range_position_temp_responses.pdf")

## Supplemental figure: p-values

ggplot(model_fits_position_fig, aes(x = p.value, fill = dir)) +
  geom_histogram(bins = 15) + 
  facet_grid(range_direction ~ term) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), labels = c(0.001, 0.01, 0.1, 1)) +
  labs(x = "Effect p-value",
       y = "Count",
       fill = "") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black")) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("dodgerblue3", "firebrick2"))

ggsave("figures/area_sensitivity/range_position_temp_pvals.pdf")

#### Strongest species response for each predictor ####

# MODO - deltaED

dove <-  clim_hab_poptrend_z %>%
  filter(Common_name == "Mourning dove")

dove_plot <- ggplot(dove, aes(x = deltaED, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Change in edge density", y = "Abundance trend", title = "Mourning dove")

# EAWP - deltaProp

pewee <- clim_hab_poptrend_z %>%
  filter(Common_name == "Eastern wood-pewee") 

pewee_plot <- ggplot(pewee, aes(x = deltaProp, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Change in forest cover", y = "", title = "Eastern wood-pewee")
 
# EATO - tmax

towhee <- clim_hab_poptrend_z %>%
  filter(Common_name == "Eastern towhee")

towhee_plot <- ggplot(towhee, aes(x = tmax, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Trend in Tmax", y = "", title = "Eastern towhee")

# RBWO - tmin

woodpecker <- clim_hab_poptrend_z %>%
  filter(Common_name == "Red-bellied woodpecker")

woodpecker_plot <- ggplot(woodpecker, aes(x = tmin, y = abundTrend)) + geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Trend in Tmin", y = "Abundance trend", title = "Red-bellied woodpecker")

# SUTA - deltaED:tmax

tanager <- clim_hab_poptrend_z %>%
  filter(Common_name == "Summer tanager") 

tanager$tmax_sign <- ifelse(tanager$tmax < 0, "-2.39 < Z-Tmax < 0", "0 < Z-Tmax < 1.57")
tanager <- ggplot(tanager, aes(x = deltaED, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~tmax_sign) +
  theme(panel.spacing = unit(4, "lines")) +
  labs(x = "Change in edge density", y = "Abundance trend", title = "Summer tanager")

# CACH - deltaED:tmin

chickadee <- clim_hab_poptrend_z %>%
  filter(Common_name == "Carolina chickadee") 

chickadee$tmin_sign <- ifelse(chickadee$tmin < 0, "-1.93 < Z-Tmin < 0", "0 < Z-Tmin < 2.30")
chickadee <- ggplot(chickadee, aes(x = deltaED, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~tmin_sign) +
  theme(panel.spacing = unit(4, "lines")) +
  labs(x = "Change in edge density", y = "Abundance trend", title = "Carolina chickadee")

## Figure for MS

indiv_spp_add <- plot_grid(dove_plot, pewee_plot, woodpecker_plot, towhee_plot, nrow = 2, labels = c("A", "B", "C", "D"))
indiv_spp_multi <- plot_grid(indiv_spp_add, tanager, chickadee, nrow = 3, rel_heights = c(0.5, 0.25, 0.25),
                             labels = c(" ", "E", "F"))
ggsave("figures/area_sensitivity/indiv_spp_multipanel.pdf", indiv_spp_multi, units = "in", width = 8, height = 12)


#### Species by site matrix ####
spp_site <- clim_hab_poptrend_mixedmod %>%
  dplyr::select(stateroute, abundTrend, SPEC) %>%
  spread(SPEC, abundTrend)
spp_cor <- cor(spp_site, use = "pairwise.complete.obs")

pdf("figures/area_sensitivity/species_correlations.pdf", height = 12, width = 12)
corrplot(spp_cor, method = "circle", tl.col = "black")
dev.off()

#### Joint models ####
obs_size <- abund_trend %>% 
  group_by(SPEC) %>% 
  summarize(nRoutes = n_distinct(stateroute)) %>%
  arrange(nRoutes) # Smallest = CERW @ 38 routes, next smallest = CAWA @ 48, NOWA @ 54

clim_hab_poptrend_mixedmod <- clim_hab_poptrend_z %>%
  filter(SPEC != "CERW")
#write.csv(clim_hab_poptrend_mixedmod, "\\\\BioArk\\hurlbertlab\\DiCecco\\data\\clim_hab_poptrend_mixedmod.csv", row.names = F)

# Mixed models
clim_hab_poptrend_mixedmod <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\clim_hab_poptrend_mixedmod.csv", stringsAsFactors = F)
clim_hab_poptrend_mixedmod <- read.csv("/Volumes/hurlbertlab/DiCecco/data/clim_hab_poptrend_mixedmod.csv", stringsAsFactors = F)

randomslope_add <- lme(abundTrend ~ tmin + tmax + ppt + deltaED + deltaProp + Wintering_Slimit_general, 
                       random = (~deltaProp + deltaED|SPEC), data = clim_hab_poptrend_mixedmod)

randomslope_fullinter <- lme(abundTrend ~ ppt + deltaProp + tmin*deltaED + tmax*deltaED + Wintering_Slimit_general, 
                         random = (~tmin + deltaED|SPEC), data = clim_hab_poptrend_mixedmod)

randomslope_noppt <- lme(abundTrend ~ tmin*deltaED + tmax*deltaED + Wintering_Slimit_general,
                               random = (~deltaProp + deltaED|SPEC), data = clim_hab_poptrend_mixedmod)

aic_table <- data.frame(model = c("Additive", "Full interactive", "Temperature and edge density"),
                        AIC = c(AIC(randomslope_add), AIC(randomslope_fullinter), AIC(randomslope_noppt)))

# Make plots of effects

add_table <- summary(randomslope_add)$tTable
add_df <- as.data.frame(add_table)
add_df$Variable <- row.names(add_table)

full_table <- summary(randomslope_fullinter)$tTable
full_df <- as.data.frame(full_table)
full_df$Variable <- row.names(full_table)

noppt_table <- summary(randomslope_noppt)$tTable
noppt_df <- as.data.frame(noppt_table)
noppt_df$Variable <- row.names(noppt_table)

model <- data.frame(model = c(rep("Additive", nrow(add_table)), 
                                    rep("Full interactive", nrow(full_table)),
                                    rep("Temperature and edge density", nrow(noppt_table))))

fixed_effects <- bind_rows(add_df, full_df, noppt_df) %>%
  bind_cols(model) %>%
  filter(Variable != "(Intercept)", !grepl("Wintering", Variable))

ggplot(fixed_effects, aes(x = Variable, y = Value)) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = Value - 1.96*Std.Error, ymax = Value + 1.96*Std.Error)) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  facet_wrap(~model) +
  theme_bw()
ggsave("figures/area_sensitivity/mixed_mods_fixedeffects.pdf", units = "in", height = 6, width = 16)

wintering_effects <- bind_rows(add_df, full_df, noppt_df) %>%
  bind_cols(model) %>%
  filter(grepl("Wintering", Variable)) %>%
  filter(model == "Full interactive")
wintering_effects$label <- c("Central America", "Mexico", "South America", "US")

ggplot(wintering_effects, aes(x = label, y = Value)) +
  geom_point(cex = 4) + geom_errorbar(aes(ymin = Value - 1.96*Std.Error, ymax = Value + 1.96*Std.Error), 
                                      width = 0, cex = 2) +
  geom_hline(yintercept = 0, lty = 2, col = "red", cex = 2) +
  theme(axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 25)) +
  theme(axis.line = element_line(colour = 'black', size = 2)) +
  labs(x = "Southern limit of wintering grounds", y = "Estimate")
ggsave("figures/area_sensitivity/mixed_mods_wintering.tiff")


# Population level responses with species responses for tmin, deltaED

clim_hab_poptrend_mixedmod$popmean <- predict(randomslope_fullinter, level = 0)
clim_hab_poptrend_mixedmod$sppmean <- predict(randomslope_fullinter, level = 1)

# Marginal plots

# Tmin
ggplot(clim_hab_poptrend_mixedmod, aes(x = tmin, color = SPEC)) + 
  geom_point(aes(y = sppmean, color = SPEC), alpha = 0.1) +
  geom_line(aes(y = sppmean, color = SPEC), alpha = 0.5, stat="smooth", method = "lm", se = F) + 
  geom_smooth(aes(y = popmean), color = "black", cex = 1, method = "lm", se = F) +
  scale_color_viridis_d() +
  labs(x = "Trend in Tmin", y = "Abundance trend", color = "Species") +
  theme(legend.position = "none")
ggsave("figures/area_sensitivity/mixedmod_randomeffect_tmin.tiff")
ggsave("figures/area_sensitivity/mixedmod_randomeffect_tmin.pdf")

# deltaED
ggplot(clim_hab_poptrend_mixedmod, aes(x = deltaED)) + 
  geom_point(aes(y = sppmean), color = "gray", alpha = 0.75) +
  geom_line(aes(y = sppmean, group = SPEC), alpha = 0.75, color = "steelblue", stat="smooth", method = "lm", se = F, cex = 0.75) + 
  geom_smooth(aes(y = popmean), color = "black", cex = 1.2, method = "lm", se = F) +
  theme(axis.line = element_line(colour = 'black', size = 1)) +
  labs(x = "Change in forest edge density", y = "Abundance trend", color = "Species") +
  theme(legend.position = "none", text = element_text(size = 24))
ggsave("figures/area_sensitivity/mixedmod_randomeffect_deltaED.tiff")
ggsave("figures/area_sensitivity/mixedmod_randomeffect_deltaED.pdf")

# Interaction of Tmax and deltaED
clim_hab_poptrend_mixedmod$dEDsign <- ifelse(clim_hab_poptrend_mixedmod$deltaED > 0, "Increase in forest edge density", "Decrease in forest edge density")

ggplot(clim_hab_poptrend_mixedmod, aes(x = tmax, color = SPEC)) + 
  geom_point(aes(y = sppmean), color = "gray", alpha = 0.75) +
  geom_line(aes(y = sppmean, group = SPEC), alpha = 0.75, color = "steelblue", stat="smooth", method = "lm", se = F, cex = 0.75) + 
  geom_smooth(aes(y = popmean), color = "black", cex = 1.2, method = "lm", se = F) +
  theme(axis.line = element_line(colour = 'black', size = 1)) +
  labs(x = "Trend in Tmax", y = "Abundance trend", color = "Species") +
  facet_wrap(~dEDsign) +
  theme(legend.position = "none") +
  theme(legend.position = "none", text = element_text(size = 24))
ggsave("figures/area_sensitivity/mixedmod_randomeffect_interaction.tiff", units = "in", height = 6, width = 12)
ggsave("figures/area_sensitivity/mixedmod_randomeffect_interaction.pdf", units = "in", height = 6, width = 12)

## Trait analysis
# Predict slope of correlation of spp abundance trends for tmin, tmax, deltaED with three traits
# Mean forest cover of routes observed on
# Mean forest edge density of routes observed on
# Climatic niche breadth

ctrl <- lmeControl(opt='optim')
randomslope_fullinter <- lme(abundTrend ~ ppt + deltaProp + tmin*deltaED + tmax*deltaED + Wintering_Slimit_general, 
                             random = (~tmin*deltaED + tmax*deltaED|SPEC), data = clim_hab_poptrend_mixedmod, control = ctrl)

z <- function(x, na.rm = TRUE) {(x - mean(x, na.rm = na.rm))/sd(x, na.rm = na.rm)}

spp_breadths_z <- data.frame(spp_breadths,
                             propFor_z = z(spp_breadths$propFor),
                             ed_z = z(spp_breadths$ed),
                             volume_z = z(spp_breadths$volume))

# Calculate relationship between abundance trend and three env predictors from model for each spp
# Join with spp env traits (z scores)
# Remake these plots!

ranef <- randomslope_fullinter$coefficients$random
ranef <- as.data.frame(ranef$SPEC)
ranef$SPEC <- row.names(ranef)

env_trends_spp <- ranef %>%
  left_join(spp_breadths_z)

# Compare random effects with indiv model estimates

compare_mods <- model_fits %>%
  dplyr::select(SPEC, term, estimate) %>%
  spread(key = "term", value = "estimate") %>%
  dplyr::select(SPEC, tmin, deltaED) %>%
  left_join(ranef, by = c("SPEC"))

ggplot(compare_mods, aes(x = tmin.x, y = tmin.y, label = SPEC)) +
  geom_text() +
  geom_abline(yintercept = 0, slope = 1) +
  labs(x = "Tmin estimate indiv. models", y = "Tmin estimate mixed model")
ggsave("figures/area_sensitivity/compare_indiv_mixedmod_tmin.pdf")

ggplot(compare_mods, aes(x = deltaED.x, y = deltaED.y, label = SPEC)) +
  geom_text() +
  geom_abline(yintercept = 0, slope = 1) +
  labs(x = "deltaED estimate indiv. models", y = "deltaED estimate mixed model")
ggsave("figures/area_sensitivity/compare_indiv_mixedmod_deltaED.pdf")


# Pairwise comparisons

neg_tmin <- filter(env_trends_spp, tmin < 0)
pos_tmin <- filter(env_trends_spp, tmin > 0)

neg_ed <- filter(env_trends_spp, deltaED < 0)
pos_ed <- filter(env_trends_spp, deltaED > 0)

shapiro.test(env_trends_spp$volume) # non-normal
shapiro.test(env_trends_spp$propFor) # non-normal
shapiro.test(env_trends_spp$ed) # normal

wilcox.test(neg_tmin$volume, pos_tmin$volume)
wilcox.test(neg_tmin$propFor, pos_tmin$propFor)
t.test(neg_tmin$ed, pos_tmin$ed)

boxplot_deltaED <- env_trends_spp %>%
  mutate(group = ifelse(deltaED < 0, "Negative", "Positive")) %>%
  left_join(breedingRange, by = c("aou" = "AOU"))

t.test(filter(boxplot_deltaED, group == "Negative")$Brange_Area_km2, 
            filter(boxplot_deltaED, group == "Positive")$Brange_Area_km2)

propFor_lm <- ggplot(boxplot_deltaED, aes(x = group, y = propFor)) + 
  geom_violin(aes(fill = group), bw = 0.1, trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1.5) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.75, cex = 2) +
  theme(legend.position = "none") + 
  ylim(c(0,1)) +
  theme(text = element_text(size = 24), axis.line = element_line(colour = 'black', size = 1),
        axis.text.x = element_text(size = 30)) +
  labs(y = "Mean proportion of forest cover", x = "Fragmentation effect") + 
  scale_x_discrete(labels = c("Negative" = "-", "Positive" = "+"))

vol_lm <- ggplot(boxplot_deltaED, aes(x = group, y = volume)) + 
  geom_violin(aes(fill = group), trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1.5) +
  scale_fill_viridis_d(begin = 0.5) +
  ylim(c(0, 4)) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.75, cex = 2) +
  theme(legend.position = "none") + 
  theme(text = element_text(size = 24), axis.line = element_line(colour = 'black', size = 1),
        axis.text.x = element_text(size = 30)) +
  labs(y = "Climatic niche breadth", x = "Fragmentation effect") + 
  scale_x_discrete(labels = c("Negative" = "-", "Positive" = "+"))

plot_grid(vol_lm, propFor_lm,
          labels = c("*", "*"),
          label_x = c(0.725, 0.74),
          label_y = c(1, 1),
          label_size = 24)
ggsave("figures/area_sensitivity/deltaED_ranef_violinplots_traits.tiff")
ggsave("figures/area_sensitivity/deltaED_ranef_violinplots_traits.pdf")

# 3 models, 1 for each slope, slope ~ for_cov + ed + volume

tmin_mod <- lm(tmin ~ propFor_z + ed_z + volume_z, data = env_trends_spp)
summary(tmin_mod)
plot(env_trends_spp$propFor, env_trends_spp$tmin)

deltaED_mod <- lm(deltaED ~ propFor_z + ed_z + volume_z, data = env_trends_spp)
plot(env_trends_spp$propFor, env_trends_spp$deltaED)

## BRMS model

### Need to use adapt_delta bigger than 0.8 for full model, see https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
library(brms)

climate_mod <- brm(abundTrend ~ ppt + tmin + tmax + Wintering_Slimit_general + (~ppt + tmin + tmax|SPEC),
                   data = clim_hab_poptrend_mixedmod,
                   inits = "0",
                   iter = 10000, 
                   cores = 4, control = list(adapt_delta = 0.97))

land_mod <- brm(abundTrend ~ deltaProp + deltaED + Wintering_Slimit_general + (~deltaProp + deltaED|SPEC),
                data = clim_hab_poptrend_mixedmod,
                inits = "0",
                iter = 10000, 
                cores = 4, control = list(adapt_delta = 0.97))

inter_mod <- brm(abundTrend ~ ppt + deltaProp + tmin*deltaED + tmax*deltaED + Wintering_Slimit_general + (~ppt + deltaProp + tmin*deltaED + tmax*deltaED|SPEC),
                 data = clim_hab_poptrend_mixedmod,
                 inits = "0",
                 iter = 10000, 
                 cores = 4, control = list(adapt_delta = 0.97))
save(climate_mod, land_mod, inter_mod, file = "model/brms_mods.Rdata")

brm_plot <- marginal_effects(inter_mod, re_formula = NULL)
plot(brm_plot, points = T)

summary <- summary(inter_mod)

load("/Volumes/hurlbertlab/DiCecco/data/brms_mods.Rdata")

#### Simulation for modeling approaches ####
# Simulated data set of five species

sim_data <- data.frame(routeno = c(1:1000))

sim_data$tmin <- rnorm(nrow(sim_data), mean = 0.2, sd = 0.6)
sim_data$deltaED <- rnorm(nrow(sim_data), mean = 0.3, sd = 0.8)

coefs <- data.frame(spec = c(1:5), 
                    tmin = c(0.5, 0.4, -0.7, -0.2, 0.3),
                    deltaED = c(-0.2, 0.7, -0.4, 0.5, 0.1),
                    int = c(0.4, 0.2, -0.1, -0.3, 0.5))
range <- 500
sim_trends <- coefs %>%
  group_by(spec) %>%
  nest() %>%
  mutate(sim = map(data, ~{
    params <- .
    routeno <- sample_n(sim_data, size = range, replace = F) %>%
      mutate(abund_trend = params$tmin*tmin + params$deltaED*deltaED + params$int*deltaED*tmin +
               rnorm(range, mean = 0, sd = 0.5))
  })) %>%
  dplyr::select(-data) %>%
  unnest()

# Individual species models

indiv_mod <- sim_trends %>%
  group_by(spec) %>%
  nest() %>%
  mutate(mod = map(data, ~{
    tidy(lm(abund_trend ~ tmin*deltaED, data = .))
  })) %>%
  dplyr::select(-data) %>%
  unnest() %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(spec, term, estimate) %>%
  spread(key = "term", value = "estimate") %>%
  mutate(model = "indiv") %>%
  rename(int = `tmin:deltaED`)

# Fixed effects model

fixed_mod <- lm(abund_trend ~ tmin*deltaED*factor(spec), data = sim_trends)
fixed_mod_df <- data.frame(term = row.names(data.frame(fixed_mod$coefficients)), est = fixed_mod$coefficients)
fixed_mod_1 <- data.frame(spec = 1, 
                          tmin = fixed_mod_df$est[2], 
                          deltaED = fixed_mod_df$est[3],
                          int = fixed_mod_df$est[8])
fixed_mod_est <- fixed_mod_df %>%
  slice(9:nrow(fixed_mod_df)) %>%
  mutate(spec = rep(c(2:5), 3)) %>%
  mutate(param = c(rep("tmin", 4), rep("deltaED", 4), rep("int", 4))) %>%
  dplyr::select(-term) %>%
  spread(key= param, value = est) %>%
  mutate(tmin = fixed_mod_1$tmin + tmin,
         deltaED = fixed_mod_1$deltaED + deltaED,
         int = fixed_mod_1$int + int) %>%
  bind_rows(fixed_mod_1) %>%
  mutate(model = "fixed")

# Mixed effects model

ctrl <- lmeControl(opt='optim')
mixed_mod <- lme(abund_trend ~ tmin*deltaED, random = (~tmin*deltaED|spec), 
                 data = sim_trends, control = ctrl)

mixed_mod_ranef <- data.frame(mixed_mod$coefficients$random[[1]])
mixed_mod_ranef$spec <- row.names(mixed_mod_ranef)
mixed_mod_fixed <- data.frame(mixed_mod$coefficients$fixed)

mixed_mod_est <- mixed_mod_ranef %>%
  mutate(tmin = tmin + mixed_mod_fixed$mixed_mod.coefficients.fixed[2],
         deltaED = deltaED + mixed_mod_fixed$mixed_mod.coefficients.fixed[3], 
         int = tmin.deltaED + mixed_mod_fixed$mixed_mod.coefficients.fixed[4]) %>%
  dplyr::select(spec, tmin, deltaED, int) %>%
  mutate(model = "mixed") %>%
  mutate(spec = as.numeric(spec))

## Compare methods

sim_results <- coefs %>%
  mutate(model = "true") %>%
  bind_rows(mixed_mod_est, indiv_mod, fixed_mod_est)

pdf("figures/area_sensitivity/simulated_data_models_equalRange.pdf")
par(mfrow = c(2,2))
  plot(filter(sim_results, model == "true")$tmin, 
       filter(sim_results, model == "mixed")$tmin, col = "red", cex = 1, pch = 19, 
       ylim= c(-0.7, 0.7),
       xlab = c("Estimated Tmin"), ylab = "True Tmin")
  abline(a = 0, b = 1)
  points(filter(sim_results, model == "true")$tmin,
         filter(sim_results, model == "indiv")$tmin, col = "blue", cex = 1, pch = 19)
  points(filter(sim_results, model == "true")$tmin,
         filter(sim_results, model == "fixed")$tmin, col = "green", cex = 1, pch = 19)
  legend(x = "topleft", legend = c("Mixed model", "Individual models", "Fixed model"), 
         col = c("red", "blue", "green"), pch = 19)

plot(filter(sim_results, model == "true")$deltaED, 
       filter(sim_results, model == "mixed")$deltaED, col = "red", cex = 1, pch = 19, 
       ylim= c(-0.7, 0.7),
       xlab = c("Estimated deltaED"), ylab = "True deltaED")
abline(a = 0, b = 1)
points(filter(sim_results, model == "true")$deltaED,
         filter(sim_results, model == "indiv")$deltaED, col = "blue", cex = 1, pch = 19)
points(filter(sim_results, model == "true")$deltaED,
         filter(sim_results, model == "fixed")$deltaED, col = "green", cex = 1, pch = 19)

plot(filter(sim_results, model == "true")$int, 
     filter(sim_results, model == "mixed")$int, col = "red", cex = 1, pch = 19, 
     ylim= c(-0.7, 0.7),
     xlab = c("Estimated deltaED:Tmin"), ylab = "True deltaED:Tmin")
abline(a = 0, b = 1)
points(filter(sim_results, model == "true")$int,
       filter(sim_results, model == "indiv")$int, col = "blue", cex = 1, pch = 19)
points(filter(sim_results, model == "true")$int,
       filter(sim_results, model == "fixed")$int, col = "green", cex = 1, pch = 19)

dev.off()
