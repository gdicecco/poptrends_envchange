# Build model of habitat fragmentation + climate ~ population trend

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
counts <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_counts_20170712.csv")
species <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_species_20170712.csv")

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

#setwd("/Volumes/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
#setwd("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
#us.proj <- readOGR("BCRs_contiguous_us.shp")

#us_subs_transf <- spTransform(us_subs, crs(us.proj))

#tm_shape() + tm_borders(us.proj) + tm_lines(us_subs_transf)

#plot(us.proj, col = "gray73", border = "gray73")
#plot(us_subs_transf, add = T)

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
  merge(filter(RT1.routes, stateroute %in% routes.short$stateroute), by = c("stateroute", "year")) %>%
  filter(year >= 1990, year < 2017)
# 1513 routes

#abund_trend <- counts.subs %>%
#  group_by(aou, stateroute) %>%
#  nest() %>%
#  mutate(lmFit = map(data, ~{
#    df <- .
#    df.short <- df %>%
#      dplyr::select(year, speciestotal) %>%
#      unique()
#    lm(speciestotal ~ year, df.short)
#  })) %>%
#  mutate(nObs = map_dbl(data, ~{
#    df <- .
#    nrow(df)
#  })) %>%
#  mutate(lm_broom = map(lmFit, tidy)) %>%
#  mutate(abundTrend = map_dbl(lm_broom, ~{
#    df <- .
#    df$estimate[2]
#  })) %>%
#  mutate(trendInt = map_dbl(lm_broom, ~{
#    df <- .
#    df$estimate[1]
#  })) %>%
#  mutate(trendPval = map_dbl(lm_broom, ~{
#    df <- .
#    df$p.value[2]
#  })) 

#hist(abund_trend$abundTrend)
#hist(abund_trend$nObs)

#abund_trend <- abund_trend %>%
#  filter(nObs > 9) %>%
#  dplyr::select(-data, -lmFit, -lm_broom)

#setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
#write.csv(abund_trend, "BBS_abundance_trends.csv", row.names = F)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
abund_trend <- read.csv("BBS_abundance_trends.csv", stringsAsFactors = F)

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

# Proportion of landscape deltas
# long to wide
  
route_ed <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge)/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2011` - `1992`,
            deltaED9201 = `2001` - `1992`,
            deltaED0692 = `2006` - `1992`) %>%
  mutate(pctT1 = abs(deltaED9201)/abs(deltaED)*100,
         pctT2 = abs(deltaED0692)/abs(deltaED)*100) %>%
  filter(pctT1 > 50 | pctT2 > 50) %>% # 2289 routes
  filter(stateroute %in% routes.short$stateroute) # 1497 routes

# Trait data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/")
traits <- read.csv("traits/spp_traits.csv", stringsAsFactors = F)

traits.short <- traits %>%
  dplyr::select(Common_name, aou, nHabitats1, nHabitats2, volume)

#ggplot(traits.short, aes(x = nHabitats1, y = volume)) + geom_point() + geom_smooth(method = "lm")
#ggsave("nhabitats_volume.pdf", units = "in")
#summary(lm(traits.short$volume ~ traits.short$nHabitats1))

# master data table
#clim_hab_poptrend <- abund_trend %>%
#  left_join(traits.short) %>%
#  left_join(route_ed, by = "stateroute") %>%
#  left_join(climate_wide, by = "stateroute") %>%
#  group_by(aou) %>%
#  mutate(abundTrend_z = (abundTrend - mean(abundTrend, na.rm = T))/sd(abundTrend, na.rm = T))

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
#write.csv(clim_hab_poptrend, "climate_fragmentation_traits_by_species.csv", row.names = F)

clim_hab_poptrend <- read.csv("climate_fragmentation_traits_by_species.csv", stringsAsFactors = F)

############ Build species-specific models #########
  
## Start with 10 most abundant species
abund_spp <- counts.subs %>% 
  group_by(aou) %>%
  summarize(spptotal = sum(speciestotal)) %>%
  arrange(desc(spptotal)) %>%
  slice(1:20) %>%
  left_join(species)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
ebird_hab <- read.csv("ebird_habitat_association.csv", stringsAsFactors = F)

spp_hab <- ebird_hab %>%
  group_by(Species) %>%
  count(Habitat) %>%
  arrange(Species, desc(n)) %>%
  group_by(Species) %>%
  filter(n == max(n)) %>%
  left_join(traits.short, by = c("Species" = "Common_name")) %>%
  left_join(data.frame(Habitat = c("Grasslands", "Croplands", "Urban and built up", "Deciduous broadleaf forest", "Mixed forest"),
                       legend = c("Grasslands", "Agricultural", "Urban", "Forest", "Forest")))

## Model selection with eBird data
# 4 models for each species: abund_trend_z ~ maxtemp + deltaED + maxtemp:deltaED
# AIC, effect sizes for each (glance, tidy functions)

## This requires deltaED for different class types
spp_models_indv <- clim_hab_poptrend %>%
  group_by(aou) %>%
  nest() %>%
  right_join(spp_hab) %>%
  mutate(data.subs = map(data, ~{
    df <- .
    df.short <- df %>%
      filter(legend == spp_hab$legend[spp_hab$Species == unique(df$Common_name)][1]) %>%
      mutate(dEDz = (deltaED - mean(deltaED))/sd(deltaED)) %>%
      dplyr::select(stateroute, abundTrend_z, legend, dEDz, tmax, tmin, ppt)
    df.short
  })) %>%
  mutate(dredged = map(data.subs, ~{
    df <- .
    mod <- lm(abundTrend_z ~ tmax + tmin + ppt + dEDz, df, na.action = na.fail)
    dobj <- dredge(mod)
    topmods <- get.models(dobj, subset = delta < 4)
    modlavg <- model.avg(topmods)
    summary(modlavg)})) %>%
  mutate(coefs = map(dredged, ~{
    sum <- . 
    sum$coefmat.full
  })) %>%
  mutate(importance = map(dredged, ~{
    sum <- .
    sum$importance
  }))

model_coefs_indv <- spp_models_indv %>%
  dplyr::select(aou, coefs) %>%
  mutate(coefs.df = map(coefs, ~{
    mat <- .
    vars <- rownames(mat)
    df <- data.frame(param = vars, mat)
  })) %>%
  dplyr::select(aou, coefs.df) %>%
  unnest()

model_importance_indv <- spp_models_indv %>%
  dplyr::select(aou, importance) %>%
  mutate(import.df = map(importance, ~{
    import <- .
    df <- data.frame(import)
    data.frame(param = rownames(df), importance = df$import)
  })) %>%
  dplyr::select(aou, import.df) %>%
  unnest()

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/model/")
write.csv(model_coefs_indv, "weighted_model_coefficients_subsetspecies.csv", row.names = F)
write.csv(model_importance_indv, "weighted_model_coefficient_importance_subsetspecies.csv", row.names = F)

model_coefs_indv <- read.csv("weighted_model_coefficients_subsetspecies.csv", stringsAsFactors = F)
model_importance_indv <- read.csv("weighted_model_coefficient_importance_subsetspecies.csv", stringsAsFactors = F)

## General model selection

abundant_spp <- abund_trend %>% 
  group_by(aou) %>% 
  summarize(nRoutes = n()) %>%
  filter(nRoutes > 40)
# 195 spp

route_ed_z <- route_ed %>%
  mutate(dEDz = (deltaED - mean(deltaED))/sd(deltaED))

route_env <- abund_trend %>%
  filter(aou %in% abundant_spp$aou) %>%
  left_join(climate_wide) %>%
  left_join(route_ed_z) %>%
  na.omit() %>% # remove species/routes with no data
  group_by(aou) %>%
  mutate(abundTrend_z = (abundTrend - mean(abundTrend, na.rm = T))/sd(abundTrend, na.rm = T))
# 198 spp, 1483 routes

# For each species:
## Make four models - abund_trend ~ dED, abund_trend ~ climate trends, abund_trend ~ dED + climate trends, abund_trend ~ dED:climate trends
## For each keep AIC, R^2, effect sizes & p vals

## Model averaging
model_avgs <- route_env %>%
  group_by(aou) %>%
  nest() %>%
  mutate(dredged = map(data, ~{
    df <- .
    mod <- lm(abundTrend_z ~ tmax + tmin + ppt + dEDz, df, na.action = na.fail)
    dobj <- dredge(mod)
    topmods <- get.models(dobj, subset = delta < 5)
    modlavg <- model.avg(topmods)
    summary(modlavg)})) %>%
  mutate(coefs = map(dredged, ~{
    sum <- . 
    sum$coefmat.full
  })) %>%
  mutate(importance = map(dredged, ~{
    sum <- .
    sum$importance
  }))

model_coefs <- model_avgs %>%
  dplyr::select(aou, coefs) %>%
  mutate(coefs.df = map(coefs, ~{
    mat <- .
    vars <- rownames(mat)
    df <- data.frame(param = vars, mat)
  })) %>%
  dplyr::select(aou, coefs.df) %>%
  unnest()

model_importance <- model_avgs %>%
  dplyr::select(aou, importance) %>%
  mutate(import.df = map(importance, ~{
    import <- .
    df <- data.frame(import)
    data.frame(param = rownames(df), importance = df$import)
  })) %>%
  dplyr::select(aou, import.df) %>%
  unnest()

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/model/")
write.csv(model_coefs, "weighted_model_coefficients.csv", row.names = F)
write.csv(model_importance, "weighted_model_coefficient_importance.csv", row.names = F)

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/model/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/model/")
model_coefs <- read.csv("weighted_model_coefficients.csv", stringsAsFactors = F)
model_importance <- read.csv("weighted_model_coefficient_importance.csv", stringsAsFactors = F)

##### Simple models
model_fits <- route_env %>%
  group_by(aou) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
  df <- .
  lm(abundTrend_z ~ tmax + tmin + ppt + dEDz, df, na.action = na.fail)
  })) %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    nrow(df)
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  dplyr::select(aou, nObs, lm_broom) %>%
  unnest()

species_estimates <- model_fits %>%
  group_by(aou) %>%
  filter(term != "(Intercept)", p.value < 0.05) %>%
  count()

env_estimates <- model_fits %>%
  filter(term != "(Intercept)", p.value < 0.05) %>%
  group_by(term) %>%
  count()

### Number of species responding positively/negatively to different drivers (simple linear models)

dir_trends_lms <- model_fits %>%
  filter(term != "(Intercept)") %>%
  mutate(significance = ifelse(p.value < 0.05, "sig", "nonsig")) %>%
  mutate(sign = ifelse(significance == "sig",ifelse(estimate < 0, "negative", "positive"), "na")) %>%
  group_by(term, significance, sign) %>%
  count()

theme_set(theme_classic())
ggplot(dir_trends_lms, aes(x = term, y = n, fill = sign)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Parameter", y = "Number of species") +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12, face = "bold")) +
  theme(axis.title.y = element_text(size = 12, face = "bold")) +
  scale_x_discrete(labels = c("dEDz" = "Edge density", "ppt" = "Ppt", "tmax" = "Tmax", "tmin" = "Tmin")) +
  scale_fill_manual(name = c("Estimate sign"),labels = c("na" = "Zero", "negative" = "Negative", "positive" = "Positive"), values = c("gray", "#0072B2", "#CC79A7"))

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/species_models_figs/")
ggsave("allspp_parameter_estimates_simpleLM.pdf", width = 6, height = 4)
ggsave("allspp_parameter_estimates_simpleLM.tiff", width = 6, height = 4, units = "in")

## Number of drivers species responding to (simple linear models)

n_trends_lm <- model_fits %>%
  filter(term!= "(Intercept)") %>%
  mutate(significance = ifelse(p.value < 0.05, 1, 0)) %>%
  group_by(aou) %>%
  summarize(ndrivers = sum(significance)) %>%
  group_by(ndrivers) %>%
  count()

ggplot(n_trends, aes(x = ndrivers, y = n, fill = ndrivers)) +
  geom_bar(stat = "identity") + 
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12, face = "bold")) +
  theme(axis.title.y = element_text(size = 12, face = "bold")) +
  labs(x = "Number of strong predictors", y = "Number of species")
ggsave("allspp_number_predictors_simpleLM.pdf", width = 5, height = 4)
ggsave("allspp_number_predictors_simpleLM.tiff", width = 5, height = 4, units = "in")

mixedresponse <- c(4950, 7610, 2890,7310,5930) # species with opposite abundance responses to inc tmax and tmin
spp_warming <- model_fits %>%
  filter(aou %in% species_estimates$aou, p.value < 0.05, term == "tmax" | term == "tmin", 
         !(aou %in% mixedresponse)) %>%
  mutate(abundDir = ifelse(estimate > 0, "increasing", "decreasing")) # 32 decreasing, 21 increasing

# Species that respond well to increases in warming temps, simple linear models

setwd("//BioArk/HurlbertLab/DiCecco/Data/")
correlates <- read.csv("Master_RO_Correlates_20110610.csv", stringsAsFactors = F)

null_traits <- traits.short %>%
  filter(aou %in% unique(route_env$aou), !(aou %in% spp_warming$aou), !(aou %in% mixedresponse)) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>% # 141 species
  dplyr::select(aou, Common_name, nHabitats1, volume, habitat, envniche) %>%
  left_join(dplyr::select(correlates, AOU, Brange_Area_km2), by = c("aou" = "AOU")) %>%
  rename(english_common_name = Common_name) %>%
  mutate(group = "null")

spp_warming_pos <- spp_warming %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(dplyr::select(correlates, AOU, Brange_Area_km2), by = c("aou" = "AOU")) %>%
  filter(term == "tmax" | term == "tmin", estimate > 0) %>%
  distinct(aou, Brange_Area_km2, english_common_name, nHabitats1, volume, habitat, envniche) %>%
  mutate(group = "warming_positive")

spp_warming_neg <- spp_warming %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(dplyr::select(correlates, AOU, Brange_Area_km2), by = c("aou" = "AOU")) %>%
  filter(term == "tmax" | term == "tmin", estimate < 0) %>%
  distinct(aou, english_common_name, nHabitats1, volume, habitat, envniche, Brange_Area_km2) %>%
  mutate(group = "warming_negative")

warming_traits <- bind_rows(null_traits, spp_warming_pos, spp_warming_neg)

warming_traits$group <- as.factor(warming_traits$group)

shapiro.test(warming_traits$nHabitats1)

wilcox.test(spp_warming_pos$nHabitats1, null_traits$nHabitats1)
wilcox.test(spp_warming_pos$volume, null_traits$volume)
wilcox.test(spp_warming_pos$Brange_Area_km2, null_traits$Brange_Area_km2) # not 


kruskal.test(nHabitats1 ~ group, data = warming_traits, na.action = na.omit)

boxplot_lm <- warming_traits %>%
  filter(group == "null" | group == "warming_positive")

nhab_lm <- ggplot(boxplot_lm, aes(x = group, y = nHabitats1)) + 
  geom_violin(aes(fill = group), bw = 0.75, trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Number of habitats") + 
  scale_x_discrete(labels = c("null" = "All other species", "warming_positive" = "Species increasing with Temp"))

vol_lm <- ggplot(boxplot_lm, aes(x = group, y = volume)) + 
  geom_violin(aes(fill = group), trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Environmental niche breadth") +
  scale_x_discrete(labels = c("null" = "All other species", "warming_positive" = "Species increasing with Temp"))

plot_grid(nhab_lm, vol_lm, 
          labels = c("*", "*"),
          label_x = 0.72,
          label_y = c(1, 1),
          label_size = 24)

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/species_models_figs/")
ggsave("allspp_warming_simpleLM_box.pdf")
ggsave("allspp_warming_simple_LM_box.tiff", units = "in")

############ Species-specific plots/results for model averaging ##############

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
spp_codes <- read.csv("four_letter_codes_birdspp.csv", stringsAsFactors = F)

## 9 species, specific habitat edge density

coefs_indv_sig <- model_coefs_indv %>%
  filter(Pr...z.. < 0.05, param != "(Intercept)") %>%
  mutate(confint = Std..Error*1.96) %>%
  left_join(spp_hab) %>%
  left_join(spp_codes, by = c("Species" = "COMMONNAME"))

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/")
theme_set(theme_classic())
ggplot(coefs_indv_sig, aes(x = SPEC, y = Estimate, color = param)) + 
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = Estimate - confint, ymax = Estimate + confint), width = 0.2, position = position_dodge(0.5), cex = 1) +
  scale_color_viridis_d(labels = c("Change in edge density", "Trend in Tmax", "Trend in Tmin")) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  theme(legend.position = c(0.15, 0.9), legend.title = element_blank()) +
  xlab("")
ggsave("spp_params.pdf")

importance_indv <- model_importance_indv %>%
  left_join(spp_hab) %>%
  left_join(spp_codes, by = c("Species" = "COMMONNAME"))

dedz <- ggplot(filter(importance_indv, param == "dEDz"), aes(x = fct_reorder(SPEC, importance), y = importance, color = param)) + 
  geom_point(size = 3) + 
  scale_color_viridis_d() +
  theme(legend.position = "none") + 
  labs(x = "", y = "Relative importance")
ppt <- ggplot(filter(importance_indv, param == "ppt"), aes(x = fct_reorder(SPEC, importance), y = importance, color = param)) + 
  geom_point(size = 3) + 
  scale_color_viridis_d(begin = 0.25) +
  theme(legend.position = "none") + 
  labs(x = "", y = "Relative importance")
tmax <- ggplot(filter(importance_indv, param == "tmax"), aes(x = fct_reorder(SPEC, importance), y = importance, color = param)) + 
  geom_point(size = 3) + 
  scale_color_viridis_d(begin = 0.5) +
  theme(legend.position = "none") + 
  labs(x = "", y = "Relative importance")
tmin <- ggplot(filter(importance_indv, param == "tmin"), aes(x = fct_reorder(SPEC, importance), y = importance, color = param)) + 
  geom_point(size = 3) + 
  scale_color_viridis_d(begin = 0.75) +
  theme(legend.position = "none") + 
  labs(x = "", y = "Relative importance")
  
plot_grid(tmax, dedz, tmin, ppt, nrow = 2,
          labels = c("Tmax", "Edge density", "Tmin", "Ppt"),
          label_x = c(0.065, 0, 0.065, 0.085))
ggsave("spp_relimport.pdf", width = 9, height = 6)

#### All species, landscape edge density
model_coefs_sig <- model_coefs %>%
  filter(Pr...z.. < 0.05, param != "(Intercept)") %>%
  mutate(confint = Std..Error*1.96) %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short)

setwd("//BioArk/HurlbertLab/DiCecco/Data/")
correlates <- read.csv("Master_RO_Correlates_20110610.csv", stringsAsFactors = F)

null_traits <- traits.short %>%
  filter(aou %in% unique(route_env$aou)) %>%
  left_join(correlates, by = c("aou" = "AOU"))

# Species that respond well to increases in habitat fragmentation 
dEDz <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  filter(param == "dEDz", Estimate > 0)

shapiro.test(dEDz$Estimate)
t.test(dEDz$nHabitats1, null_traits$nHabitats1)
t.test(dEDz$volume, null_traits$volume) # Marginal
t.test(dEDz$Brange_Area_km2, null_traits$Brange_Area_km2)
t.test(dEDz$NumBiomes, null_traits$NumBiomes)

# Species that respond poorly to increases in habitat fragmentation
dEDz <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  filter(param == "dEDz", Estimate < 0)

shapiro.test(dEDz$Estimate)
t.test(dEDz$nHabitats1, null_traits$nHabitats1)
t.test(dEDz$volume, null_traits$volume)
t.test(dEDz$Brange_Area_km2, null_traits$Brange_Area_km2)
t.test(dEDz$NumBiomes, null_traits$NumBiomes)


# Species that respond well to decreases in ppt
ppt <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  filter(param == "ppt", Estimate < 0)

shapiro.test(ppt$Estimate)
wilcox.test(ppt$nHabitats1, null_traits$nHabitats1) # significant
wilcox.test(ppt$volume, null_traits$volume) # marginal
wilcox.test(ppt$Brange_Area_km2, null_traits$Brange_Area_km2)
wilcox.test(ppt$NumBiomes, null_traits$NumBiomes)

# Species that respond well to increases in ppt
ppt <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  filter(param == "ppt", Estimate > 0)

shapiro.test(ppt$Estimate)
t.test(ppt$nHabitats1, null_traits$nHabitats1)
t.test(ppt$volume, null_traits$volume)
t.test(ppt$Brange_Area_km2, null_traits$Brange_Area_km2)
t.test(ppt$NumBiomes, null_traits$NumBiomes)


# Species that respond well to increases in tmin
tmin <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  filter(param == "tmin", Estimate > 0)

shapiro.test(tmin$Estimate)

wilcox.test(tmin$nHabitats1, null_traits$nHabitats1) # significant
wilcox.test(tmin$volume, null_traits$volume) # significant
wilcox.test(tmin$Brange_Area_km2, null_traits$Brange_Area_km2)
wilcox.test(tmin$NumBiomes, null_traits$NumBiomes)

# Species that respond well to increases in tmax
tmax <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  left_join(correlates, by = c("aou" = "AOU")) %>%
  filter(param == "tmax", Estimate > 0)

shapiro.test(tmax$Estimate)

t.test(tmax$nHabitats1, null_traits$nHabitats1) # marginal
t.test(tmax$volume, null_traits$volume)
t.test(tmin$Brange_Area_km2, null_traits$Brange_Area_km2)
t.test(tmin$NumBiomes, null_traits$NumBiomes)

### Number of species responding positively/negatively to different drivers

dir_trends <- model_coefs %>%
  filter(param != "(Intercept)") %>%
  mutate(confint = Std..Error*1.96) %>%
  mutate(significance = ifelse(Pr...z.. < 0.05, "sig", "nonsig")) %>%
  mutate(sign = ifelse(significance == "sig",ifelse(Estimate < 0, "negative", "positive"), "na")) %>%
  group_by(param, significance, sign) %>%
  count()

theme_set(theme_classic())
ggplot(dir_trends, aes(x = param, y = n, fill = sign)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Parameter", y = "Number of species") +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12, face = "bold")) +
  theme(axis.title.y = element_text(size = 12, face = "bold")) +
  scale_x_discrete(labels = c("dEDz" = "Edge density", "ppt" = "Ppt", "tmax" = "Tmax", "tmin" = "Tmin")) +
  scale_fill_manual(name = c("Estimate sign"),labels = c("na" = "Zero", "negative" = "Negative", "positive" = "Positive"), values = c("gray", "#0072B2", "#CC79A7"))

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/species_models_figs/")
ggsave("allspp_parameter_estimates.pdf", width = 6, height = 4)
ggsave("allspp_parameter_estimates.tiff", width = 6, height = 4, units = "in")

## Number of drivers species responding to

n_trends <- model_coefs %>%
  filter(param != "(Intercept)") %>%
  mutate(significance = ifelse(Pr...z.. < 0.05, 1, 0)) %>%
  group_by(aou) %>%
  summarize(ndrivers = sum(significance)) %>%
  group_by(ndrivers) %>%
  count()

ggplot(n_trends, aes(x = ndrivers, y = n, fill = ndrivers)) +
  geom_bar(stat = "identity") + 
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12, face = "bold")) +
  theme(axis.title.y = element_text(size = 12, face = "bold")) +
  labs(x = "Number of strong predictors", y = "Number of species")
ggsave("allspp_number_predictors.pdf", width = 5, height = 4)
ggsave("allspp_number_predictors.tiff", width = 5, height = 4, units = "in")

### Species that are responding strongly to multiple sources of change
twodrivers <- model_coefs_sig %>%
  group_by(aou) %>%
  count() %>%
  filter(n > 1)
# 14 species

twodrivers_coefs <- model_coefs_sig %>%
  filter(aou %in% twodrivers$aou)

### Box plots for species tolerant of warming in tmax and tmin

# Species that respond well to increases in tmin
tmin <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  filter(param == "tmin", Estimate > 0)

shapiro.test(tmin$Estimate)

wilcox.test(tmin$nHabitats1, null_traits$nHabitats1) # significant
wilcox.test(tmin$volume, null_traits$volume) # significant

boxplot <- null_traits %>%
  mutate(group = ifelse(null_traits$aou %in% tmin$aou, "tolerant", "allspp"))

nhab <- ggplot(boxplot, aes(x = group, y = nHabitats1)) + 
  geom_violin(aes(fill = group), trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Number of habitats") + 
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species increasing with Tmin"))
  
vol <- ggplot(boxplot, aes(x = group, y = volume)) + 
  geom_violin(aes(fill = group), trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Environmental niche width") +
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species increasing with Tmin"))
  

plot_grid(nhab, vol, 
          labels = c("*", "*"),
          label_x = 0.725,
          label_y = c(1, 1))

ggsave("allspp_tmin_box.pdf")
ggsave("allspp_tmin_box.tiff", units = "in")

# Species that respond well to increases in tmax
tmax <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  filter(param == "tmax", Estimate > 0)

shapiro.test(tmax$Estimate)

t.test(tmax$nHabitats1, null_traits$nHabitats1) # marginal
t.test(tmax$volume, null_traits$volume)

boxplot2 <- null_traits %>%
  mutate(group = ifelse(null_traits$aou %in% tmax$aou, "tolerant", "allspp"))

nhab2 <- ggplot(boxplot2, aes(x = group, y = nHabitats1)) + 
  geom_violin(aes(fill = group), trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Number of habitats") + 
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species increasing with Tmax"))

vol2 <- ggplot(boxplot2, aes(x = group, y = volume)) + 
  geom_violin(aes(fill = group), trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Environmental niche width") +
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species increasing with Tmax"))


plot_grid(nhab2, vol2)

ggsave("allspp_tmax_box.pdf")
ggsave("allspp_tmax_box.tiff", units = "in")

# Species that respond well to increases in warming temps
temp <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  filter(param == "tmax" | param == "tmin", Estimate > 0)

shapiro.test(temp$Estimate)

wilcox.test(temp$nHabitats1, null_traits$nHabitats1) # significant
wilcox.test(temp$volume, null_traits$volume) # significant

boxplot2 <- null_traits %>%
  mutate(group = ifelse(null_traits$aou %in% temp$aou, "tolerant", "allspp"))

nhab_temp <- ggplot(boxplot2, aes(x = group, y = nHabitats1)) + 
  geom_violin(aes(fill = group), bw = 0.75, trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Number of habitats") + 
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species increasing with Temp"))

vol_temp <- ggplot(boxplot2, aes(x = group, y = volume)) + 
  geom_violin(aes(fill = group), trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Environmental niche width") +
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species increasing with Temp"))

plot_grid(nhab_temp, vol_temp, 
          labels = c("*", "*"),
          label_x = 0.732,
          label_y = c(1, 1))

ggsave("allspp_warming_box.pdf")
ggsave("allspp_warming_box.tiff", units = "in")

##### Community-level responses #####

traits.short
clim_hab_poptrend
climate_trends

route_frag <- frags %>%
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge)/sum(total.area)) %>%
  spread(key = "year", value = "ED")

# Are there community differences between fragmented and unfragmented routes?
## Group routes as fragmented or intact

landEDQ <- quantile(route_ed$deltaED, c(0.33, 0.66, 1))

routes_noEDchange <- route_ed %>%
  filter(deltaED > landEDQ[1] & deltaED < landEDQ[2])

routes_frag_noChange <- route_frag %>%
  filter(stateroute %in% routes_noEDchange$stateroute) # 500 routes

## Plot: Trait distribution comparisons between routes in fragmented and intact landscapes - continuous

medianHab <- median(traits.short$nHabitats2, na.rm = T)
medianVol <- median(traits.short$volume, na.rm =T)

route_traits_cont <- clim_hab_poptrend %>%
  filter(stateroute %in% routes_noEDchange$stateroute) %>%
  left_join(routes_frag_noChange) %>%
  group_by(stateroute) %>%
  summarize(landED = mean(`2011`, na.rm = T),
            sppRich = n(),
            nHabSpec = sum(nHabitats2 < 4, na.rm = T),
            nHabGen = sum(nHabitats2 > 4, na.rm = T),
            nVolSpec = sum(volume < medianVol, na.rm = T),
            nVolGen = sum(volume > medianVol, na.rm = T),
            pHabSpec = nHabSpec/sppRich,
            pHabGen = nHabGen/sppRich,
            pVolSpec = nVolSpec/sppRich,
            pVolGen = nVolGen/sppRich)

theme_set(theme_classic())

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/community_comparisons/")

ggplot(route_traits_cont, aes(x = landED, y = sppRich)) + geom_point() + geom_smooth(method = "lm")
summary(lm(route_traits_cont$sppRich ~ route_traits_cont$landED))
ggsave("sppRich_fragmentation_continuous.pdf", units = "in")

habGrp <- ggplot(route_traits_cont, aes(x = landED, y = nHabSpec)) + geom_point() + geom_smooth(method = "lm") + ylab("Number of habitat specialists")
habGrp2 <- ggplot(route_traits_cont, aes(x = landED, y = nHabGen)) + geom_point() + geom_smooth(method = "lm") + ylab("Number of habitat generalists")
volGrp <- ggplot(route_traits_cont, aes(x = landED, y = nVolSpec)) + geom_point() + geom_smooth(method = "lm") + ylab("Number of environmental niche specialists")
volGrp2 <- ggplot(route_traits_cont, aes(x = landED, y = nVolGen)) + geom_point() + geom_smooth(method = "lm") + ylab("Number of environmental niche generalists")

plot_grid(habGrp, habGrp2, volGrp, volGrp2, nrow = 2)
ggsave("number_specialists_generalists_fragmentation.pdf", height = 8, width = 10, units = "in")

summary(lm(route_traits_cont$nHabSpec ~ route_traits_cont$landED))
summary(lm(route_traits_cont$nHabGen ~ route_traits_cont$landED))
summary(lm(route_traits_cont$nVolSpec ~ route_traits_cont$landED))
summary(lm(route_traits_cont$nVolGen ~ route_traits_cont$landED))

habGrpP <- ggplot(route_traits_cont, aes(x = landED, y = pHabSpec)) + geom_point() + geom_smooth(method = "lm") + ylab("Proportion habitat specialists") + xlab("Landscape edge density")
habGrp2P <- ggplot(route_traits_cont, aes(x = landED, y = pHabGen)) + geom_point() + ylab("Proportion habitat generalists") + xlab("Landscape edge density")
volGrpP <- ggplot(route_traits_cont, aes(x = landED, y = pVolSpec)) + geom_point() + ylab("Proportion environmental niche specialists") + xlab("Landscape edge density")
volGrp2P <- ggplot(route_traits_cont, aes(x = landED, y = pVolGen)) + geom_point() + ylab("Proportion environmental niche generalists") + xlab("Landscape edge density")

plot_grid(habGrpP, habGrp2P, volGrpP, volGrp2P, nrow = 2)
ggsave("prop_specialists_generalists_fragmentation.pdf", height = 8, width = 10, units = "in")

summary(lm(route_traits_cont$pHabSpec ~ route_traits_cont$landED))
summary(lm(route_traits_cont$pHabGen ~ route_traits_cont$landED))
summary(lm(route_traits_cont$pVolSpec ~ route_traits_cont$landED))
summary(lm(route_traits_cont$pVolGen ~ route_traits_cont$landED))


## Plot: map of these routes and their distribution in US

routelist_frag_coords <- routes_noEDchange %>%
  left_join(routes)

setwd("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/")
us_sf <- read_sf("BCRs_contiguous_us.shp")

routes_sf <- st_as_sf(routelist_frag_coords, coords = c("longitude", "latitude"))

us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "gray")

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/community_comparisons/")
route_map <- us + tm_shape(routes_sf) + 
  tm_dots(alpha = 0.5, size = 1)
route_map
tmap_save(route_map, "routes_fragmentation_noChange_map.pdf", units = "in")

## Group species into specialists/generalists
medianHab
medianVol

spp_trends_habgrps <- clim_hab_poptrend %>%
  filter(stateroute %in% routes_noEDchange$stateroute) %>%
  left_join(route_frag) %>%
  mutate(habitat = ifelse(nHabitats2 > medianHab, "generalist", "specialist"),
         envniche = ifelse(volume > medianVol, "generalist", "specialist")) %>%
  mutate(abundDir = ifelse(trendPval > 0.05, "stable", ifelse(abundTrend > 0, "increasing", "decreasing"))) %>%
  group_by(stateroute, habitat) %>%
  mutate(nSpp = n()) %>%
  group_by(stateroute, habitat, abundDir) %>%
  summarize(n = n(),
            pctN = n/unique(nSpp),
            landED = mean(`2011`, na.rm = T))

ggplot(filter(spp_trends_habgrps, !is.na(abundDir), !is.na(habitat)), aes(x = landED, y = pctN, col = abundDir)) + 
  geom_point() + geom_smooth(method = "lm", se = F) + facet_grid(~habitat) + labs(y = "Proportion of species")
ggsave("habitat_groups_abundTrends_nochangeFrag_CC.pdf", units = "in", height = 6, width = 12)

spp_trends_volgrps <- clim_hab_poptrend %>%
  filter(stateroute %in% routes_noEDchange$stateroute) %>%
  left_join(route_frag) %>%
  mutate(habitat = ifelse(nHabitats2 > medianHab, "generalist", "specialist"),
         envniche = ifelse(volume > medianVol, "generalist", "specialist")) %>%
  mutate(abundDir = ifelse(trendPval > 0.05, "stable", ifelse(abundTrend > 0, "increasing", "decreasing"))) %>%
  group_by(stateroute, envniche) %>%
  mutate(nSpp = n()) %>%
  group_by(stateroute, envniche, abundDir) %>%
  summarize(n = n(),
            pctN = n/unique(nSpp),
            landED = mean(`2011`, na.rm = T))

ggplot(filter(spp_trends_volgrps, !is.na(abundDir), !is.na(envniche)), aes(x = landED, y = pctN, col = abundDir)) + 
  geom_point() + geom_smooth(method = "lm", se = F) + facet_grid(~envniche) + labs(y = "Proportion of species")
ggsave("envniche_groups_abundTrends_nochangeFrag_CC.pdf", units = "in", height = 6, width = 12)

summary(lm(filter(spp_trends_habgrps, habitat == "generalist", abundDir == "decreasing")$pctN ~ filter(spp_trends_habgrps, habitat == "generalist", abundDir == "decreasing")$landED))
summary(lm(filter(spp_trends_habgrps, habitat == "generalist", abundDir == "increasing")$pctN ~ filter(spp_trends_habgrps, habitat == "generalist", abundDir == "increasing")$landED))
summary(lm(filter(spp_trends_habgrps, habitat == "generalist", abundDir == "stable")$pctN ~ filter(spp_trends_habgrps, habitat == "generalist", abundDir == "stable")$landED))

summary(lm(filter(spp_trends_habgrps, habitat == "specialist", abundDir == "decreasing")$pctN ~ filter(spp_trends_habgrps, habitat == "specialist", abundDir == "decreasing")$landED))
summary(lm(filter(spp_trends_habgrps, habitat == "specialist", abundDir == "increasing")$pctN ~ filter(spp_trends_habgrps, habitat == "specialist", abundDir == "increasing")$landED))
summary(lm(filter(spp_trends_habgrps, habitat == "specialist", abundDir == "stable")$pctN ~ filter(spp_trends_habgrps, habitat == "specialist", abundDir == "stable")$landED))

summary(lm(filter(spp_trends_volgrps, envniche == "generalist", abundDir == "decreasing")$pctN ~ filter(spp_trends_volgrps, envniche == "generalist", abundDir == "decreasing")$landED))
summary(lm(filter(spp_trends_volgrps, envniche == "generalist", abundDir == "increasing")$pctN ~ filter(spp_trends_volgrps, envniche == "generalist", abundDir == "increasing")$landED))
summary(lm(filter(spp_trends_volgrps, envniche == "generalist", abundDir == "stable")$pctN ~ filter(spp_trends_volgrps, envniche == "generalist", abundDir == "stable")$landED))

summary(lm(filter(spp_trends_volgrps, envniche == "specialist", abundDir == "decreasing")$pctN ~ filter(spp_trends_volgrps, envniche == "specialist", abundDir == "decreasing")$landED))
summary(lm(filter(spp_trends_volgrps, envniche == "specialist", abundDir == "increasing")$pctN ~ filter(spp_trends_volgrps, envniche == "specialist", abundDir == "increasing")$landED))
summary(lm(filter(spp_trends_volgrps, envniche == "specialist", abundDir == "stable")$pctN ~ filter(spp_trends_volgrps, envniche == "specialist", abundDir == "stable")$landED))



# Routes with no change in fragmentation but climate change - 289 routes total

## Routes with significant climate change

clim_sig_trends <- climate_trends %>%
  filter(trendPval < 0.05) %>%
  distinct(stateroute)

routes_climchange_n <- clim_hab_poptrend %>%
  filter(stateroute %in% clim_sig_trends$stateroute) %>%
  filter(stateroute %in% routelist_frag$stateroute) %>%
  left_join(route_frag) %>%
  mutate(abundDir = ifelse(trendPval > 0.05, "stable", ifelse(abundTrend > 0, "increasing", "decreasing"))) %>%
  group_by(stateroute, abundDir) %>%
  summarize(n = n(),
            landED = mean(`2011`, na.rm = T))

## Plot: Abundance trends at routes with no change in fragmentation but climate change

ggplot(filter(routes_climchange_n, !is.na(abundDir)), aes(x = landED, y = n, col = abundDir)) + geom_point() + geom_smooth(method = "lm")
ggsave("abund_trend_inc_dec_fragmentation.pdf", units = "in")

## Map of these routes

routes_climchange_nochangeF <- clim_hab_poptrend %>%
  filter(stateroute %in% clim_sig_trends$stateroute) %>%
  filter(stateroute %in% routelist_frag$stateroute) %>%
  left_join(route_frag) %>%
  left_join(routes) %>%
  dplyr::select(stateroute, latitude, longitude, deltaED, ppt, tmax, tmin, `2011`, -aou)

hist(routes_climchange_nochangeF$`2011`)
hist(routes_climchange_nochangeF$tmin)
hist(routes_climchange_nochangeF$ppt)
hist(routes_climchange_nochangeF$tmax)

routes_CC_sf <- st_as_sf(routes_climchange_nochangeF, coords = c("longitude", "latitude")) %>%
  dplyr::rename(landED = `2011`)

us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "gray")

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/community_comparisons/")
route_map <- us + tm_shape(routes_CC_sf) + 
  tm_dots(size = 0.5, col = "landED")
route_map
tmap_save(route_map, "routes_nochangeFragmentation_ClimChange_map.pdf", units = "in")

# Change in trait space (direction) on routes with increases in fragmentation

## Routes with increase in fragmentation

inc_frag_routes <- route_ed %>%
  filter(deltaED > landEDQ[2]) %>%
  left_join(routes)

hist(inc_frag_sf$deltaED)

## Map of these routes:

inc_frag_sf <- st_as_sf(inc_frag_routes, coords = c("longitude", "latitude"))

us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "gray")

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/community_comparisons/")
route_map <- us + tm_shape(inc_frag_sf) + 
  tm_dots(size = 0.5, alpha = 0.5)
route_map
tmap_save(route_map, "routes_inc_fragmentation_map.pdf", units = "in")

## Need communities in 1990-1995 and communities in 2011-2016

traits.short$habitat <- ifelse(traits.short$nHabitats2 > medianHab, "generalist", "specialist")
traits.short$envniche <- ifelse(traits.short$volume > medianVol, "generalist", "specialist")

route_communities <- counts.subs %>%
  filter(stateroute %in% inc_frag_routes$stateroute) %>%
  filter(year %in% c(1990:1995, 2011:2016)) %>%
  mutate(time_window = ifelse(year < 1996, "t1", "t2")) %>%
  group_by(stateroute, time_window) %>%
  distinct(aou) %>%
  left_join(traits.short) %>%
  left_join(dplyr::select(inc_frag_routes, stateroute, deltaED, latitude, longitude))

bothT <- route_communities %>% 
  group_by(stateroute) %>% 
  distinct(time_window) %>% 
  count() %>% 
  filter(n == 2)

route_communities_both <- route_communities %>%
  filter(stateroute %in% bothT$stateroute)

## Change in % generalists on routes with increases in fragmentation

route_percent_change <- route_communities_both %>%
  group_by(stateroute, deltaED, time_window) %>%
  summarize(sppRich = n(),
            pctHabGen = sum(habitat == "generalist", na.rm = T)/n(),
            pctVolGen = sum(envniche == "generalist", na.rm = T)/n()) %>%
  group_by(stateroute, deltaED) %>%
  summarize(dSpRich = sppRich[time_window == "t2"] - sppRich[time_window == "t1"],
            dPctHab = pctHabGen[time_window == "t2"] - pctHabGen[time_window == "t1"],
            dPctVol = pctVolGen[time_window == "t2"] - pctVolGen[time_window == "t1"]) %>%
  left_join(climate_wide)

pcthab <- ggplot(route_percent_change, aes(x = deltaED, y = dPctHab)) + geom_point() + geom_smooth(method = "lm", se = F) +
  ylab("Change in % habitat generalists")
pctvol <- ggplot(route_percent_change, aes(x = deltaED, y = dPctVol)) + geom_point() + geom_smooth(method = "lm", se = F) +
  ylab("Change in % environmental niche generalists")
plot_grid(pcthab, pctvol, nrow = 1)
ggsave("deltaPercentGeneralists_inc_fragmentation.pdf", height = 6, width = 14, units = "in")

ggplot(route_percent_change, aes(x = deltaED, y = dSpRich)) + geom_point() + geom_smooth(method = "lm", se = F) +
  ylab("Change in species richness")
ggsave("inc_fragmentation_deltaSppRichness.pdf", units = "in")

sppRichMod <- glm(dSpRich ~ deltaED + ppt + tmax + tmin + deltaED:tmax + deltaED:tmin, data = route_percent_change)
habGenMod <- summary(lm(dPctHab ~ deltaED, data = route_percent_change))
volGenMod <- summary(lm(dPctVol ~ deltaED, data = route_percent_change)) # dec. in deltaED predicts inc. in dPctVol

####### Community-level responses: raster heat maps ########

## Stable fragmentation amt routes

# Plot: number of habitats vs. edge density, color = % of species in each category
route_traits_habbins <- clim_hab_poptrend %>%
  filter(stateroute %in% routes_noEDchange$stateroute, !is.na(nHabitats2)) %>%
  left_join(routes_frag_noChange) %>%
  group_by(stateroute) %>%
  mutate(sppRich = n()) %>%
  group_by(stateroute, nHabitats2) %>%
  summarize(landED = mean(`2011`, na.rm = T),
            propS = n()/mean(sppRich))

ggplot(route_traits_habbins, aes(x = landED, y = factor(nHabitats2), z = propS)) +
  stat_summary_2d(bins = 15) + scale_fill_viridis_c() +
  labs(fill = "Prop. of species",
       x = "Landscape edge density",
       y = "Number of breeding habitats")
ggsave("figures/community_comparisons/heatplot_nhab_propS.pdf", units = "in", height = 8, width = 10)

# Plot: number of habitats vs. edge density, color = slope of abundance trend
boxplot(clim_hab_poptrend$abundTrend[clim_hab_poptrend$abundTrend > -30])

route_traits_hab_abund <- clim_hab_poptrend %>%
  filter(stateroute %in% routes_noEDchange$stateroute, !is.na(nHabitats2)) %>%
  filter(abundTrend > -30) %>%
  left_join(routes_frag_noChange) %>%
  group_by(stateroute, nHabitats2) %>%
  summarize(landED = mean(`2011`, na.rm = T),
            meanAT = mean(abundTrend))

ggplot(route_traits_hab_abund, aes(x = landED, y = factor(nHabitats2), z = meanAT)) +
  stat_summary_2d(bins = 15) + scale_fill_viridis_c() +
  labs(fill = "Mean abund. trend",
       x = "Landscape edge density",
       y = "Number of breeding habitats")
ggsave("figures/community_comparisons/heatplot_nhab_abundTrend.pdf", units = "in", height = 8, width = 10)

# Plot: env niche breadth vs. edge density, color = % of species in each category
route_traits_vol <- clim_hab_poptrend %>%
  filter(stateroute %in% routes_noEDchange$stateroute) %>%
  left_join(routes_frag_noChange) %>%
  filter(!is.na(volume)) %>%
  mutate(vol_bin = 0.2*floor(volume/0.2) + 0.2/2) %>%
  group_by(stateroute) %>%
  mutate(sppRich = n()) %>%
  group_by(stateroute, vol_bin) %>%
  summarize(landED = mean(`2011`, na.rm = T),
            propS = n()/mean(sppRich))

ggplot(route_traits_vol, aes(x = landED, y = factor(vol_bin), z = propS)) +
    stat_summary_2d(bins = 15) + scale_fill_viridis_c() +
    labs(fill = "Prop. of species",
         x = "Landscape edge density",
         y = "Environmental niche breadth")
ggsave("figures/community_comparisons/heatplot_vol_propS.pdf", units = "in", height = 8, width = 10)
  
  
# Plot: env niche breadth vs. edge density, color = slope of abundance trend

route_traits_vol_abund <- clim_hab_poptrend %>%
  filter(stateroute %in% routes_noEDchange$stateroute) %>%
  filter(!is.na(volume)) %>%
  filter(abundTrend > -30) %>%
  mutate(vol_bin = 0.2*floor(volume/0.2) + 0.2/2) %>%
  left_join(routes_frag_noChange) %>%
  group_by(stateroute, vol_bin) %>%
  summarize(landED = mean(`2011`, na.rm = T),
            meanAT = mean(abundTrend))

ggplot(route_traits_vol_abund, aes(x = landED, y = factor(vol_bin), z = meanAT)) +
  stat_summary_2d(bins = 15) + scale_fill_viridis_c() +
  labs(fill = "Mean abund. trend",
       x = "Landscape edge density",
       y = "Environmental niche breadth")
ggsave("figures/community_comparisons/heatplot_vol_abundTrend.pdf", units = "in", height = 8, width = 10)

##### Community-level responses: forest fragmentation #####

## Filter out only species who use forests

habitats <- read.csv("traits/spp_habitat_detailed.csv", stringsAsFactors = F)

forest_spp <- habitats %>%
  group_by(aou) %>%
  filter(habitat1 == "Forest") %>%
  distinct(aou) # 273 species that use forest

forest_traits <- traits.short %>%
  filter(aou %in% forest_spp$aou) %>% # 234 spp
  filter(Common_name != "unknown") # 219 spp

counts.forest <- counts.subs %>%
  filter(aou %in% forest_traits$aou)

poptrend_forest <- clim_hab_poptrend %>%
  filter(aou %in% forest_traits$aou)

# Forest fragmentation at routes

forest_ed <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge[legend == "Forest"])/sum(total.area)) %>%
  spread(key = "year", value = "ED")

forest_deltaED <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge[legend == "Forest"])/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2011` - `1992`,
            deltaED9201 = `2001` - `1992`,
            deltaED0692 = `2006` - `1992`) %>%
  mutate(pctT1 = abs(deltaED9201)/abs(deltaED)*100,
         pctT2 = abs(deltaED0692)/abs(deltaED)*100) %>%
  filter(pctT1 > 50 | pctT2 > 50) %>%
  filter(stateroute %in% routes.short$stateroute) # 1483 routes

## ID routes with no change in forest fragmentation 

forestEDQ <- quantile(forest_deltaED$deltaED, c(0.33, 0.66, 1))

forest_noEDchange <- forest_deltaED %>%
  filter(deltaED > forestEDQ[1] & deltaED < forestEDQ[2])

forest_frag_noChange <- forest_ed %>%
  filter(stateroute %in% forest_noEDchange$stateroute) %>% # 489 routes
  rename("forestED" = `2011`) %>%
  dplyr::select(stateroute, forestED)

## Community composition at routes with no change in forest fragmentation 

forest_communities <- counts.forest %>%
  left_join(forest_traits) %>%
  filter(stateroute %in% forest_frag_noChange$stateroute) %>%
  left_join(forest_frag_noChange) %>%
  group_by(stateroute) %>%
  distinct(aou, nHabitats2, volume, forestED) 

# Plot: number of breeding habitats vs. forest edge density, color = proportion of species in each habitat bin
forest_habitat_bins <- forest_communities %>%
  group_by(stateroute) %>%
  mutate(sppRich = n()) %>%
  group_by(stateroute, nHabitats2) %>%
  summarize(forestED = mean(forestED, na.rm = T),
            propS = n()/mean(sppRich))

ggplot(filter(forest_habitat_bins, !is.na(nHabitats2)), aes(x = forestED, y = factor(nHabitats2), z = propS)) +
  stat_summary_2d(bins = 15) + scale_fill_viridis_c() +
  labs(fill = "Prop. of species",
       x = "Forest edge density",
       y = "Number of breeding habitats")
ggsave("figures/community_comparisons_forestED/heatplot_nhab_propS.pdf", units = "in", height = 8, width = 10)

# Plot: mean number of breeding habitats of route

forest_habitat_mean <- forest_communities %>%
  group_by(forestED) %>%
  summarize(meanHB = mean(nHabitats2, na.rm = T),
            stdevHB = sd(nHabitats2, na.rm = T))

ggplot(forest_habitat_mean, aes(x = forestED, y = meanHB)) + 
  geom_errorbar(aes(ymin = meanHB - stdevHB, ymax = meanHB + stdevHB), alpha = 0.4) +
  geom_point() + geom_smooth(method = "loess") +
  labs(x = "Forest edge density", y = "Mean number of breeding habitats")
ggsave("figures/community_comparisons_forestED/forestED_vs_meanhabitatbreadth.pdf", units = "in")

# Things to do:
## Consider cropping observations geographically (eastern deciduous forests)
## Create scale bar for forest edge density routes
## Try analysis for proportion urban as well?
