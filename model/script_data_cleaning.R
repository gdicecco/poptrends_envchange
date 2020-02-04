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

### Read in data #####

# Bird population data
## BBS 2017 Version
routes <- read.csv("\\\\BioArk\\HurlbertLab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\BioArk\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# species four letter codes
spp_codes <- read.csv("traits/four_letter_codes_birdspp.csv", stringsAsFactors = F)
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
# write.csv(fourletter_codes, "traits/four_letter_codes_aou.csv", row.names = F)

## Population trends
counts.subs <- counts %>%
  filter(aou %in% area_aous$aou) %>%
  merge(filter(RT1.routes, stateroute %in% routes.short$stateroute), by = c("stateroute", "year")) %>%
  filter(rpid == 101) %>%
  filter(year >= 1990, year < 2017)
# 1834 routes

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
# 1990-2010 for Canada, 1992-2016 for US
abund_trend <- counts.subs %>%
  left_join(obs_years, by = c('stateroute', 'year')) %>%
  group_by(aou, stateroute) %>%
  nest() %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    country <- unique(df$countrynum.x)
    if(country == 840){
      df <- df %>%
        filter(year >= 1992 & year <= 2016) %>%
        unique()
      nrow(df)
    } else {
      df <- df %>%
        filter(year >= 1990 & year <= 2010) %>%
        unique()
      nrow(df)
    }
  })) %>%
  filter(nObs > 9) %>%
  mutate(lmFit = purrr::map(data, ~{
    df <- .
    country <- unique(df$countrynum.x)
    if(country == 840){
      df.short <- df %>%
        dplyr::select(year, first_yr, speciestotal) %>%
        filter(year >= 1992 & year <= 2016) %>%
        unique()
      glm(speciestotal ~ year + first_yr, family = poisson, data = df.short)
    } else {
      df.short <- df %>%
        dplyr::select(year, first_yr, speciestotal) %>%
        filter(year >= 1990 & year <= 2010) %>%
        unique()
      glm(speciestotal ~ year + first_yr, family = poisson, data = df.short)
    }
  }))  %>%
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
  dplyr::select(-data, -lmFit, -lm_broom)

# write.csv(abund_trend, "model/BBS_abundance_trends.csv", row.names = F)

abund_trend <- read.csv("model/BBS_abundance_trends.csv", stringsAsFactors = F) %>%
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
temp_range <- read.csv("traits/bbs_aou_temp_range.csv", stringsAsFactors = F)

setwd("//BioArk/HurlbertLab/DiCecco/Data/")
correlates <- read.csv("Master_RO_Correlates_20110610.csv", stringsAsFactors = F)

breedingRange <- correlates %>%
  dplyr::select(AOU, Brange_Area_km2)

# Climate data
climate_trends <- read.csv("climate/bbs_routes_climate_trends.csv", stringsAsFactors = F)

climate_wide <- climate_trends %>%
  dplyr::select(-trendPval) %>%
  spread(key = "env", value = "climateTrend")

# Habitat fragmentation data
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
  left_join(env_change, by = "stateroute") %>%
  filter(!is.na(propForest))

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

