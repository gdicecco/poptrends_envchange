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
area_sensitive <- read.csv("traits/area-sensitivity-forest-birds-Bouliner1998.csv", stringsAsFactors = F)
area_sensitive <- area_sensitive[, -c(4)]

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

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
abund_trend <- read.csv("BBS_abundance_trends.csv", stringsAsFactors = F) %>%
  filter(aou %in% area_aous$aou) %>% ## Only ones lost from area_aous are nocturnal: owls, nightjars (Caprimulgus)
  left_join(area_aous) %>%
  filter(sporder != "Accipitriformes") ## Remove hawks

# Trait data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")
setwd("/Users/gracedicecco/Desktop/git/NLCD_fragmentation/")
traits <- read.csv("traits/spp_traits.csv", stringsAsFactors = F)

traits.short <- traits %>%
  dplyr::select(Common_name, aou, nHabitats1, nHabitats2, volume)

setwd("//BioArk/HurlbertLab/DiCecco/Data/")
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
  filter(year == 2011, propForest > 0.25) %>% # use filter to ID routes
  arrange(ED)

forest_deltaED <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  group_by(stateroute, year) %>%
  summarize(ED = sum(total.edge[legend == "Forest"])/sum(total.area)) %>%
  spread(key = "year", value = "ED") %>%
  group_by(stateroute) %>%
  summarize(deltaED = `2011` - `1992`)

forest_deltaP <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(propForest = prop.landscape) %>%
  spread(key = "year", value = "propForest") %>%
  group_by(stateroute) %>%
  summarize(deltaProp = `2011` - `1992`)

forest <- forest_ed %>%
  left_join(forest_deltaED) %>%
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

#### Route trends ####

# Plot of routes used

forest_routes <- us_subs[us_subs$rteno %in% forest_ed$stateroute, ]

us.proj <- readOGR("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/BCRs_contiguous_us.shp")

us_subs_transf <- spTransform(forest_routes, crs(us.proj))

plot(us.proj, col = "gray73", border = "gray73")
plot(us_subs_transf, add = T)


hist(clim_hab_poptrend$deltaED)
hist(clim_hab_poptrend$deltaProp)

cor(clim_hab_poptrend$ED, clim_hab_poptrend$propForest) # -0.41
cor(clim_hab_poptrend$deltaED, clim_hab_poptrend$deltaProp) # -0.25
cor(dplyr::select(clim_hab_poptrend, ppt, tmin, tmax, deltaED, deltaProp, ED, propForest))

# Breadth of forest ED and forest cover for area and non-area sensitive species

env_breadth <- clim_hab_poptrend %>%
  dplyr::group_by(Area_sensitivity, SPEC, aou) %>%
  summarize(ed = mean(ED),
            std_ed = sd(ED),
            propFor = mean(propForest),
            std_for = sd(propForest),
            patchArea = mean(meanPatchArea),
            std_patch = sd(meanPatchArea))

ed_breadth <- ggplot(env_breadth, aes(x = reorder(SPEC, ed), y = ed, color = SPEC)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ed - 1.96*std_ed, ymax = ed + 1.96*std_ed)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Forest edge density") +
  scale_color_viridis_d() +
  theme(legend.position = "none")

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
  filter(year == 2011)

forest_deltaP_allroutes <- frags %>% # 2314 routes
  left_join(newcode, by = c("class" = "code")) %>%
  filter(legend == "Forest") %>%
  group_by(stateroute, year) %>%
  summarize(propForest = prop.landscape) %>%
  spread(key = "year", value = "propForest") %>%
  group_by(stateroute) %>%
  summarize(deltaProp = `2011` - `1992`)

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

volume <- traits %>%
  filter(aou %in% env_breadth_allroutes$aou) %>%
  left_join(dplyr::select(env_breadth_allroutes, SPEC, aou)) %>%
  dplyr::select(SPEC, aou, volume)

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

spp_breadths <- env_breadth_allroutes %>%
  dplyr::select(SPEC, aou, propFor, min_for, max_for) %>%
  left_join(dplyr::select(env_breadth, SPEC, aou, ed, std_ed)) %>%
  left_join(dplyr::select(traits, aou, volume))

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

cor(spp_breadths[, c(4,7,9)]) # edxpropFor = -0.65, edxVolume = 0.11, volumexPropFor = -0.59

#### Individual models ####
# of climate + frag + loss + climate:frag + climate:loss for species

z <- function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T)}

forest_z <- data.frame(stateroute = forest$stateroute,
                       ED = z(forest$ED),
                       propForest = z(forest$propForest),
                       deltaED = z(forest$deltaED), 
                       deltaProp = z(forest$deltaProp))

climate_wide_z <- climate_wide %>%
  mutate_at(.vars = c("ppt", "tmin", "tmax"),
            .funs = z)

clim_hab_poptrend_z <- abund_trend %>%
  left_join(dplyr::select(traits.short, -Common_name)) %>%
  left_join(forest_z, by = "stateroute") %>%
  left_join(climate_wide_z, by = "stateroute") %>%
  filter(!is.na(propForest)) %>%
  group_by(SPEC, aou) %>% 
  mutate(abundTrend_z = (abundTrend - mean(abundTrend, na.rm = T))/sd(abundTrend, na.rm = T)) %>%
  left_join(fourletter_codes)

# Models by species

model_fits <- clim_hab_poptrend_z %>%
  group_by(aou, SPEC) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    lm(abundTrend_z ~ tmax*deltaED + tmin*deltaED + tmax*deltaProp + tmin*deltaProp, df, na.action = na.fail)
  })) %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    nrow(df)
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  dplyr::select(aou, SPEC, nObs, lm_broom) %>%
  unnest() %>%
  filter(nObs > 40) %>%
  filter(term != "(Intercept)")

range(model_fits$nObs) # 1 spp at 38, only one below 40

# Distribution of effect sizes for species

density_plot <- function(variable, label) {
  ggplot(filter(model_fits, term == variable), aes(x = estimate)) + 
    geom_density(fill = "honeydew3", color = "honeydew4", size = 1) + 
    geom_vline(aes(xintercept = mean(estimate)), lty = 2, cex = 1) + 
    labs(x = label, y = "Number of species")
}

deltaED <- density_plot("deltaED", "Change in ED")
deltaProp <- density_plot("deltaProp", "Change in forest cover")
tmin <- density_plot("tmin", "Trend in Tmin")
tmax <- density_plot("tmax", "Trend in Tmax")
tmaxED <- density_plot("tmax:deltaED", "Tmax:change in ED")
tminED <- density_plot("deltaED:tmin", "Tmin:change in ED")
tmaxProp <- density_plot("tmax:deltaProp", "Tmax:change in forest cover")
tminProp <- density_plot("tmin:deltaProp", "Tmin:change in forest cover")

plot_grid(deltaED, deltaProp, tmin, tmax, nrow = 2)
ggsave("figures/area_sensitivity/indiv_effects_distributions.pdf", units = "in", height = 6, width = 8)

plot_grid(tmaxED, tminED, tmaxProp, tminProp, nrow = 2)
ggsave("figures/area_sensitivity/indiv_interact_distributions.pdf", units = "in", height = 6, width = 8)

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
  mutate(trait_mod = map(data, ~{
    df <- .
    lm(estimate ~  0 + volume + propFor + Trophic.Group + Foraging + migclass, data = df)
  }),
  tidy = map(trait_mod, ~{
    mod <- .
    tidy(mod)
  }),
  r2 = map_dbl(trait_mod, ~{
    mod <- .
    glance(mod)$r.squared
  })) %>%
  dplyr::select(term, tidy, r2) %>%
  unnest()

spp_traits_pred <- spp_traits %>%
  filter(p.value < 0.1)

#### Make a joint model: fixed effects[climate + forest + winteringGround] + random slopes/intercepts[species]
# No z score abundance trends, use term for southern limit of wintering grounds

obs_size <- abund_trend %>% 
  group_by(SPEC) %>% 
  summarize(nRoutes = n_distinct(stateroute)) %>%
  arrange(nRoutes) # Smallest = CERW @ 38 routes, next smallest = CAWA @ 48, NOWA @ 54

clim_hab_poptrend_mixedmod <- clim_hab_poptrend_z %>%
  filter(SPEC != "CERW")
#write.csv(clim_hab_poptrend_mixedmod, "\\\\BioArk\\hurlbertlab\\DiCecco\\data\\clim_hab_poptrend_mixedmod.csv", row.names = F)

# Species by site matrix
spp_site <- clim_hab_poptrend_mixedmod %>%
  dplyr::select(stateroute, abundTrend, SPEC) %>%
  spread(SPEC, abundTrend)
spp_cor <- cor(spp_site, use = "pairwise.complete.obs")

pdf("figures/area_sensitivity/species_correlations.pdf", height = 12, width = 12)
corrplot(spp_cor, method = "circle", tl.col = "black")
dev.off()

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

### Simulate data to compare modeling approaches:
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
