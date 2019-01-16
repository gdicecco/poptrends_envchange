# Build model of habitat fragmentation + climate ~ population trend

library(tidyverse)
library(raster)
library(rgdal)

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

library(tmap)
tm_shape() + tm_borders(us.proj) + tm_lines(us_subs_transf)

plot(us.proj, col = "gray73", border = "gray73")
plot(us_subs_transf, add = T)

routes.short <- RT1.routes %>% # subset stateroutes that were filtered by criteria above
  filter(stateroute %in% us_subs@data$rteno)
# 2161 routes

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
  group_by(aou, stateroute) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    df.short <- df %>%
      dplyr::select(year, speciestotal) %>%
      unique()
    lm(speciestotal ~ year, df.short)
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
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    df$p.value[2]
  })) 

hist(abund_trend$abundTrend)
hist(abund_trend$nObs)

abund_trend <- abund_trend %>%
  filter(nObs > 9) %>%
  dplyr::select(-data, -lmFit, -lm_broom)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
write.csv(abund_trend, "BBS_abundance_trends.csv", row.names = F)

abund_trend <- read.csv("BBS_abundance_trends.csv", stringsAsFactors = F)

# Climate data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
climate_trends <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)

climate_wide <- climate_trends %>%
  dplyr::select(-trendPval) %>%
  spread(key = "env", value = "climateTrend")

# Habitat fragmentation data
frags <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_nlcd.csv", stringsAsFactors = F)

classlegend00s <- data.frame(class = c(11:12, 21:24, 31, 41:43, 51, 52, 71, 81:82, 90, 95), 
                             legend = c("Open water", "Perennial ice/snow", "Developed, open space", "Developed, low intensity", "Developed, medium intensity", "Developed, high intensity", "Barren land", "Deciduous forest", "Evergreen forest", "Mixed forest", "Dwarf scrub", "Shrub/scrub", "Grassland/herbaceous", "Pasture/hay", "Cultivated crops", "Woody wetlands", "Emergent herbaceous wetlands"))
classlegend92 <- data.frame(class = c(11:12, 85, 21:23, 31:33, 41:43, 51, 61, 71, 81:84, 91:92),
                            legend = c("Open water", "Perennial ice/snow", "Developed, open space", "Developed, low intensity", "Developed, medium intensity", "Developed, high intensity", "Barren land", "Barren land", "Barren land", "Deciduous forest", "Evergreen forest", "Mixed forest", "Shrub/scrub", "Cultivated crops", "Grassland/herbaceous", "Pasture/hay",  "Cultivated crops",  "Cultivated crops",  "Cultivated crops", "Woody wetlands", "Emergent herbaceous wetlands"))

frags.92 <- frags %>%
  filter(year == 1992) %>%
  left_join(classlegend92)
frags.00s <- frags %>%
  filter(year > 1992) %>%
  left_join(classlegend00s)

# Filter out land cover classes of interest only
# Proportion of landscape deltas
# long to wide
frags.legend <- bind_rows(frags.92, frags.00s) %>%
  dplyr::select(year, stateroute, edge.density, legend)  %>%
  spread(key = "year", value = "edge.density")
colnames(frags.legend) <-  c("stateroute", "legend", "ED1992", "ED2001", "ED2006", "ED2011")

frag_trends <- frags.legend %>%
  replace_na(list(ED1992 = 0, ED2001 = 0, ED2006 = 0, ED2011 = 0)) %>%
  mutate(deltaED = ED2011 - ED1992)

# Trait data
traits <- read.csv("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/traits/spp_traits.csv", stringsAsFactors = F)

############ Build models #########
  
traits.short <- traits %>%
  dplyr::select(Common_name, aou, nHabitats1, nHabitats2, volume)

# master data table
clim_hab_poptrend <- abund_trend %>%
  left_join(traits.short) %>%
  left_join(frag_trends, by = "stateroute") %>%
  left_join(climate_wide, by = "stateroute")
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
write.csv(clim_hab_poptrend, "climate_fragmentation_traits_by_species.csv", row.names = F)
clim_hab_poptrend <- read.csv("climate_fragmentation_traits_by_species.csv", stringsAsFactors = F)

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
                       legend = c("Grassland/herbaceous", "Cultivated crops", "Developed, medium intensity", "Deciduous forest", "Mixed forest")))

## Model selection with eBird data
# 4 models for each species: abund_trend ~ maxtemp + deltaED + maxtemp:deltaED
# AIC, effect sizes for each (glance, tidy functions)

library(MuMIn)
spp_models_indv <- clim_hab_poptrend %>%
  group_by(aou) %>%
  nest() %>%
  right_join(spp_hab) %>%
  mutate(data.subs = map(data, ~{
    df <- .
    df.short <- df %>%
      filter(legend == spp_hab$legend[Species == unique(df$Common_name)]) %>%
      mutate(dEDz = (deltaED - mean(deltaED))/sd(deltaED)) %>%
      dplyr::select(stateroute, abundTrend, legend, dEDz, tmax, tmin, ppt)
    df.short
  })) %>%
  mutate(dredged = map(data.subs, ~{
    df <- .
    mod <- lm(abundTrend ~ tmax + tmin + ppt + dEDz, df, na.action = na.fail)
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

route_ed <- frags %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge)/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2011` - `1992`) %>%
  mutate(dEDz = (deltaED - mean(deltaED))/sd(deltaED))

abundant_spp <- abund_trend %>% 
  group_by(aou) %>% 
  summarize(nRoutes = n()) %>%
  filter(nRoutes > 40)
# 198 spp

route_env <- abund_trend %>%
  filter(aou %in% abundant_spp$aou) %>%
  left_join(climate_wide) %>%
  left_join(route_ed) %>%
  na.omit() # remove species/routes with no data
# 198 spp, 1690 routes

# For each species:
## Make four models - abund_trend ~ dED, abund_trend ~ climate trends, abund_trend ~ dED + climate trends, abund_trend ~ dED:climate trends
## For each keep AIC, R^2, effect sizes & p vals

model_avgs <- route_env %>%
  group_by(aou) %>%
  nest() %>%
  mutate(dredged = map(data, ~{
    df <- .
    mod <- lm(abundTrend ~ tmax + tmin + ppt + dEDz, df, na.action = na.fail)
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

model_coefs <- read.csv("weighted_model_coefficients.csv", stringsAsFactors = F)
model_importance <- read.csv("weighted_model_coefficient_importance.csv", stringsAsFactors = F)

############ Plots/Results ##############

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")

spp_codes <- read.csv("four_letter_codes_birdspp.csv", stringsAsFactors = F)

## 9 species, specific habitat edge density

coefs_indv_sig <- model_coefs_indv %>%
  filter(Pr...z.. < 0.05, param != "(Intercept)") %>%
  mutate(confint = Std..Error*1.96) %>%
  left_join(spp_hab) %>%
  left_join(spp_codes, by = c("Species" = "COMMONNAME"))

library(ggplot2)
library(cowplot)

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

library(forcats)
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

## All species, landscape edge density
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
t.test(dEDz$volume, null_traits$volume)
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
t.test(ppt$nHabitats1, null_traits$nHabitats1) # marginal
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

wilcox.test(tmin$nHabitats1, null_traits$nHabitats1) 
wilcox.test(tmin$volume, null_traits$volume) # marginal
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

t.test(tmax$nHabitats1, null_traits$nHabitats1)
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

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/figures/")
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

### Species that are responding strongly to temperature change and land use change
twodrivers <- model_coefs_sig %>%
  group_by(aou) %>%
  count() %>%
  filter(n > 1)
# 11 species

twodrivers_coefs <- model_coefs_sig %>%
  filter(aou %in% twodrivers$aou)

viridis_palette <- c("#440154FF", "#33638DFF", "#3CBB75FF", "#FDE725FF")
pospos <- ggplot(filter(twodrivers_coefs, aou == 4400 | aou == 4060 | aou == 5010), aes(x = english_common_name, y = Estimate, color = param)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = Estimate - confint, ymax = Estimate + confint), width = 0.2, position = position_dodge(0.5), cex = 1) +
  scale_color_manual(values = viridis_palette, labels = c("dEDz" = "Edge density", "ppt" = "Ppt", "tmax" = "Tmax", "tmin" = "Tmin")) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  theme(legend.position = c(0.1, 0.9), legend.title = element_blank()) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 12)) +
  ylim(c(-9.5, 28))
pospos

negneg <- ggplot(filter(twodrivers_coefs, aou == 6250 | aou == 5190 | aou == 5810 | aou == 6870 | aou == 6730),
                 aes(x = english_common_name, y = Estimate, color = param)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = Estimate - confint, ymax = Estimate + confint), width = 0.2, position = position_dodge(0.5), cex = 1) +
  scale_color_manual(values = viridis_palette, labels = c("dEDz" = "Edge density", "ppt" = "Ppt", "tmax" = "Tmax", "tmin" = "Tmin")) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  theme(legend.position = c(0.1, 0.9), legend.title = element_blank()) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 12)) +
  ylim(c(-9.5, 28))
negneg

negpos <- ggplot(filter(twodrivers_coefs, aou == 7310), aes(x = english_common_name, y = Estimate, color = param)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = Estimate - confint, ymax = Estimate + confint), width = 0.2, position = position_dodge(0.5), cex = 1) +
  scale_color_manual(values = viridis_palette, labels = c("dEDz" = "Edge density", "ppt" = "Ppt", "tmax" = "Tmax", "tmin" = "Tmin")) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  theme(legend.position = c(0.1, 0.9), legend.title = element_blank()) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 12)) +
  ylim(c(-9.5, 28))
negpos

legend <- get_legend(negneg)

plot_grid(pospos + theme(legend.position = "none"),
          negneg + theme(legend.position = "none"),
          negpos + theme(legend.position = "none"),
          legend,
          nrow = 1,
          labels = c("Benefitting from change", "Vulnerable to change", "Mixed response"),
          hjust = c(-0.32, -0.2, -0.25), 
          rel_widths = c(1,1,1,0.3))
ggsave("allspp_twodrivers.pdf", width = 16, height = 4, units = "in")
ggsave("allspp_twodrivers.tiff", width = 16, height = 4, units = "in")

### Box plots for species tolerant of warming in tmax and tmin

# Species that respond well to increases in tmin
tmin <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  filter(param == "tmin", Estimate > 0)

shapiro.test(tmin$Estimate)

wilcox.test(tmin$nHabitats1, null_traits$nHabitats1) 
wilcox.test(tmin$volume, null_traits$volume) # marginal

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
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species benefitted by inc. Tmin"))
  
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
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species benefitted by inc. Tmin"))
  

plot_grid(nhab, vol, 
          labels = c("", "*"),
          label_x = 0.725,
          label_y = c(1, 0.94))

ggsave("allspp_tmin_box.pdf")
ggsave("allspp_tmin_box.tiff", units = "in")

# Species that respond well to increases in tmax
tmax <- model_coefs_sig %>%
  left_join(species) %>%
  left_join(spp_codes, by = c("english_common_name" = "COMMONNAME")) %>%
  left_join(traits.short) %>%
  filter(param == "tmax", Estimate > 0)

shapiro.test(tmax$Estimate)

t.test(tmax$nHabitats1, null_traits$nHabitats1)
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
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species benefitted by inc. Tmax"))

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
  scale_x_discrete(labels = c("allspp" = "All species", "tolerant" = "Species benefitted by inc. Tmax"))


plot_grid(nhab2, vol2)

ggsave("allspp_tmax_box.pdf")
ggsave("allspp_tmax_box.tiff", units = "in")

