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

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)
setwd("\\\\BioArk/hurlbertlab/Databases/BBS/GPS_stoplocations/")
us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data$rteno > 3000), ]

routes.short <- routes %>% # subset stateroutes that were filtered by criteria above
  filter(stateroute %in% us_subs@data$rteno)

# Subset species
landbirds <- species %>%
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010)

## Population trends
counts.subs <- counts %>%
  filter(aou %in% landbirds$aou) %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  filter(year > 1990, year < 2017)

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

## Start with 10 most abundant species
abund_spp <- counts.subs %>% 
  group_by(aou) %>%
  summarize(spptotal = sum(speciestotal)) %>%
  arrange(desc(spptotal)) %>%
  slice(1:20) %>%
  left_join(species)

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

spp_models <- clim_hab_poptrend %>%
  group_by(aou) %>%
  nest() %>%
  right_join(spp_hab) %>%
  mutate(data.subs = map(data, ~{
    df <- .
    df.short <- df %>%
      filter(legend == spp_models$legend[Species == unique(df$Common_name)]) %>%
      dplyr::select(stateroute, abundTrend, legend, deltaED, tmax)
    df.short
  })) %>%
  mutate(clim_mod = map(data.subs, ~{
    df <- .
    lm(abundTrend ~ tmax, df)
  })) %>%
  mutate(ed_mod = map(data.subs, ~{
    df <- .
    lm(abundTrend ~ deltaED, df)
  })) %>%
  mutate(add_mod = map(data.subs, ~{
    df <- .
    lm(abundTrend ~ tmax + deltaED, df)
  })) %>%
  mutate(int_mod = map(data.subs, ~{
    df <- .
    lm(abundTrend ~ tmax + deltaED + tmax*deltaED, df)
  })) %>%
  mutate(clim_tidy = map(clim_mod, tidy)) %>%
  mutate(clim_glance = map(clim_mod, glance)) %>%
  mutate(ed_tidy = map(ed_mod, tidy)) %>%
  mutate(ed_glance = map(ed_mod, glance)) %>%
  mutate(add_tidy = map(add_mod, tidy)) %>%
  mutate(add_glance = map(add_mod, glance)) %>%
  mutate(int_tidy = map(int_mod, tidy)) %>%
  mutate(int_glance = map(int_mod, glance))
  

spp_mod_glances <- spp_models %>%
  dplyr::select(aou, Species, Habitat, nHabitats1, nHabitats2, volume, legend, clim_glance, ed_glance, add_glance, int_glance) %>%
  gather(key = model, value = glance, clim_glance, ed_glance, add_glance, int_glance) %>%
  unnest() %>%
  arrange(aou, AIC)

spp_mod_params <- spp_models %>%
  dplyr::select(aou, Species, Habitat, nHabitats1, nHabitats2, volume, legend, clim_tidy, ed_tidy, add_tidy, int_tidy) %>%
  gather(key = model, value = tidy, clim_tidy, ed_tidy, add_tidy, int_tidy) %>%
  unnest() %>%
  arrange(aou, model, term)

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

route_env <- abund_trend %>%
  filter(aou %in% abundant_spp$aou) %>%
  left_join(climate_wide) %>%
  left_join(route_ed)

# For each species:
## Make four models - abund_trend ~ dED, abund_trend ~ cliamte trends, abund_trend ~ dED + climate trends, abund_trend ~ dED:climate trends
## For each keep AIC, R^2, effect sizes & p vals

models <- route_env %>%
  group_by(aou) %>%
  nest() %>%
  mutate(clim_mod = map(data, ~{
    df <- .
    lm(abundTrend ~ tmax + tmin + ppt, df)
  })) %>%
  mutate(ed_mod = map(data, ~{
    df <- .
    lm(abundTrend ~ dEDz, df)
  })) %>%
  mutate(add_mod = map(data, ~{
    df <- .
    lm(abundTrend ~ tmax + tmin + ppt + dEDz, df)
  })) %>%
  mutate(int_mod = map(data, ~{
    df <- .
    lm(abundTrend ~ tmax + tmin + ppt + dEDz + tmax*dEDz + tmin*dEDz + ppt*dEDz, df)
  })) %>%
  mutate(clim_tidy = map(clim_mod, tidy)) %>%
  mutate(clim_glance = map(clim_mod, glance)) %>%
  mutate(ed_tidy = map(ed_mod, tidy)) %>%
  mutate(ed_glance = map(ed_mod, glance)) %>%
  mutate(add_tidy = map(add_mod, tidy)) %>%
  mutate(add_glance = map(add_mod, glance)) %>%
  mutate(int_tidy = map(int_mod, tidy)) %>%
  mutate(int_glance = map(int_mod, glance))

model_glances <- models %>%
  dplyr::select(aou, clim_glance, ed_glance, add_glance, int_glance) %>%
  gather(key = model, value = glance, clim_glance, ed_glance, add_glance, int_glance) %>%
  unnest() %>%
  arrange(aou, AIC) %>%
  group_by(aou) %>%
  mutate(dAIC = AIC - min(AIC))

model_params <- models %>%
  dplyr::select(aou, clim_tidy, ed_tidy, add_tidy, int_tidy) %>%
  gather(key = model, value = tidy, clim_tidy, ed_tidy, add_tidy, int_tidy) %>%
  unnest() %>% 
  filter(term != "(Intercept)") %>%
  arrange(aou, model, term)

best_mods <- model_glances %>%
  filter(dAIC < 2)

library(MuMIn)


model_avgs <- route_env %>%
  group_by(aou) %>%
  nest() %>%
  mutate(dredged = map(data, ~{
    df <- .
    mod <- lm(abundTrend ~ tmax + tmin + ppt + dEDz + tmax*dEDz + tmin*dEDz + ppt*dEDz, df, na.action = na.fail)
    dobj <- dredge(mod)
    modlavg <- model.avg(dobj)
    summary(modlavg)})) %>%
  mutate(coefs = map(dredged, ~{
    sum <- . 
    sum$coefmat.subset
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

############ Plots ##############
