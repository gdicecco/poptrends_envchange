### Fit species-specific models using area-sensitivity classification

### Libraries

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

### Read in data

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

### Route trends

hist(clim_hab_poptrend$deltaED)
hist(clim_hab_poptrend$deltaProp)

cor(clim_hab_poptrend$ED, clim_hab_poptrend$propForest) # -0.41
cor(clim_hab_poptrend$deltaED, clim_hab_poptrend$deltaProp) # -0.25
cor(dplyr::select(clim_hab_poptrend, ppt, tmin, tmax, deltaED, deltaProp, ED, propForest))

# Are area-sensitive species found more on routes with greater proportion of forest and less fragmentation?

prop_spp <- clim_hab_poptrend %>%
  group_by(stateroute) %>%
  mutate(nSpp = n_distinct(aou)) %>%
  group_by(stateroute, propForest, ED, deltaED, deltaProp, Area_sensitivity) %>%
  summarize(propSpp = n_distinct(aou)/mean(nSpp))

stable <- prop_spp %>%
  filter(deltaProp < 0.05, deltaProp > -0.05, deltaED < 0.05, deltaED > -0.05)

edFor <- ggplot(stable, aes(x = ED, y = propSpp, color = Area_sensitivity)) + 
  geom_point() + geom_smooth(se= F) +
  scale_color_viridis_d(begin = 0.5) +
  labs(x = "Forest edge density", y = "Proportion of species", color = "Area sensitivity")

pFor <- ggplot(stable, aes(x = propForest, y = propSpp, color = Area_sensitivity)) + 
  geom_point() + geom_smooth(se= F) +
  scale_color_viridis_d(begin = 0.5) +
  theme(legend.position = "none") +
  labs(x = "Proportion forest cover", y = "Proportion of species", color = "Area sensitivity")

plot_grid(pFor, edFor, rel_widths = c(0.42, 0.58))
ggsave("figures/area_sensitivity/proportion_species.pdf", units = "in", width = 12, height = 5)

### Are there habitat/environmental niche differences between area-sensitive, non-area-sensitive species?

area_traits <- clim_hab_poptrend %>%
  distinct(aou, nHabitats1, nHabitats2, volume, Area_sensitivity)

hist(area_traits$nHabitats2)
shapiro.test(area_traits$nHabitats2)

wilcox.test(filter(area_traits, Area_sensitivity == "Area-sensitive")$nHabitats2, 
            filter(area_traits, Area_sensitivity == "Non-area-sensitive")$nHabitats2)

hist(area_traits$volume)
shapiro.test(area_traits$volume)

wilcox.test(filter(area_traits, Area_sensitivity == "Area-sensitive")$volume, 
            filter(area_traits, Area_sensitivity == "Non-area-sensitive")$volume)

hab <- ggplot(area_traits, aes(x = Area_sensitivity, y = nHabitats2, fill = Area_sensitivity)) + 
  geom_violin(trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Number of breeding habitats")

vol <- ggplot(area_traits, aes(x = Area_sensitivity, y = volume, fill = Area_sensitivity)) + 
  geom_violin(trim = F, draw_quantiles = c(0.5), alpha = 0.5, cex = 1) +
  scale_fill_viridis_d(begin = 0.5) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5) +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12)) +
  labs(y = "Environmental niche breadth")

plot_grid(hab, vol,
labels = c("*", ""),
label_x = 0.735,
label_y = c(1, 1))
ggsave("figures/area_sensitivity/trait_distribution.pdf", units = "in")

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

cover_breadth <- ggplot(env_breadth_allroutes, aes(x = reorder(SPEC, propFor), y = propFor, color = SPEC)) +
  geom_point() + 
  scale_color_viridis_d()+
  geom_errorbar(aes(ymin = min_for, ymax = max_for)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Forest cover") +
  theme(legend.position = "none")

volume <- traits %>%
  filter(aou %in% env_breadth_allroutes$aou) %>%
  left_join(dplyr::select(env_breadth_allroutes, SPEC, aou)) %>%
  dplyr::select(SPEC, aou, volume)

volume_plot <- ggplot(volume, aes(x = reorder(SPEC, volume), y = volume, color = SPEC)) +
  geom_point() + 
  scale_color_viridis_d()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Climatic niche breadth") +
  theme(legend.position = "none")

plot_grid(ed_breadth, cover_breadth, volume_plot, nrow = 2)
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
ggsave("figures/area_sensitivity/traits_covariation.pdf", units = "in", height = 8, width = 8)

cor(spp_breadths[, c(4,7,9)]) # edxpropFor = -0.65, edxVolume = 0.11, volumexPropFor = -0.59

### Individual models of climate + frag + loss + climate:frag + climate:loss for species

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
  mutate(abundTrend_z = (abundTrend - mean(abundTrend, na.rm = T))/sd(abundTrend, na.rm = T))

# Models by species

model_fits <- clim_hab_poptrend %>%
  group_by(aou, Area_sensitivity) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    lm(abundTrend_z ~ tmax + tmin + ppt + deltaED + deltaProp, df, na.action = na.fail)
  })) %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    nrow(df)
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  dplyr::select(aou, Area_sensitivity, nObs, lm_broom) %>%
  unnest()

range(model_fits$nObs) # 1 spp at 38, only one below 40

# Comparison of parameter estimates for area-sensitive vs. non-area-sensitive spp
## For each parameter, plot of confidence intervals around estimate for each species, color those by area sensitivity

for(i in 2:6) {
param <- unique(model_fits$term)[[i]]
ggplot(filter(model_fits, term == param), aes(x = reorder(aou, estimate), y = estimate, color = Area_sensitivity)) +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error)) +
    scale_color_viridis_d(begin = 0.5) +
    labs(x = "Species", y = "Estimate", title = param) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank()) +
    facet_wrap(~Area_sensitivity)
ggsave(paste0("figures/area_sensitivity/", param, "_estimates.pdf"), units = "in", height = 6, width = 10)
}

# Use z-score abundance trends, not env. variables
theme_set(theme_bw())

ggplot(clim_hab_poptrend_z, aes(x = deltaED, y = abundTrend_z, color = factor(aou))) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "lm", se = F) +
  scale_color_viridis_d() + facet_wrap(~Area_sensitivity) +
  theme(legend.position = "none") +
  labs(x = "Change in forest edge density", y = "Abundance trend")
ggsave("figures/area_sensitivity/abundance_trends_deltaED.pdf", units = "in")

ggplot(clim_hab_poptrend_z, aes(x = deltaProp, y = abundTrend_z, color = factor(aou))) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "lm", se = F) +
  scale_color_viridis_d() + facet_wrap(~Area_sensitivity) + 
  theme(legend.position = "none") +
  labs(x = "Change in proportion of forest cover", y = "Abundance trend")
ggsave("figures/area_sensitivity/abundance_trends_deltaProp.pdf", units = "in")

spp_slopes_both <- clim_hab_poptrend_z %>%
  group_by(Area_sensitivity, aou) %>%
  nest() %>%
  mutate(lm = map(data, ~{
    df <- .
    tidy(lm(abundTrend_z ~ deltaProp + deltaED, data = df))
  })) %>%
  dplyr::select(Area_sensitivity, aou, lm) %>%
  unnest() %>%
  mutate(upper = estimate + 1.96*std.error, 
         lower = estimate - 1.96*std.error)

ggplot(filter(spp_slopes_both, term == "deltaProp"), 
       aes(x = reorder(factor(aou), estimate), y = estimate, color = Area_sensitivity)) +
  geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, size = 1, lty = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "AOU", y = "Parameter estimate: abundance trend ~ change in forest cover")
ggsave("figures/area_sensitivity/abundance_slopes_deltaProp.pdf", units = "in")

ggplot(filter(spp_slopes_both, term == "deltaED"), 
       aes(x = reorder(factor(aou), estimate), y = estimate, color = Area_sensitivity)) +
  geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, size = 1, lty = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "AOU", y = "Parameter estimate: abundance trend ~ change in forest ED")
ggsave("figures/area_sensitivity/abundance_slopes_deltaED.pdf", units = "in")

# Variance partitioning of change in forest cover and edge density for each species

spp_varpart <- clim_hab_poptrend_z %>%
  group_by(Area_sensitivity, aou) %>%
  nest() %>%
  mutate(both = map(data, ~{
    df <- .
    summary(lm(abundTrend_z ~ deltaProp + deltaED, data = df))$r.squared
  })) %>%
  mutate(propFor = map(data, ~{
    df <- .
    summary(lm(abundTrend_z ~ deltaProp, data = df))$r.squared
  })) %>%
  mutate(ED = map(data, ~{
    df <- .
    summary(lm(abundTrend_z ~ deltaED, data = df))$r.squared
  })) %>%
  dplyr::select(Area_sensitivity, aou, both, propFor, ED) %>%
  unnest() %>%
  mutate(forest_cover = both - ED,
         forest_ED = both - propFor,
         shared = propFor - forest_cover) %>%
  dplyr::select(Area_sensitivity, aou, forest_cover, forest_ED, shared) %>%
  gather(key = "variable", value = "variance", c(3:5)) %>%
  group_by(aou) %>%
  mutate(total = sum(variance))

ggplot(spp_varpart, aes(x = reorder(aou, total), y = variance, fill = variable)) +
  geom_col(position = "stack") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "AOU", y = "Variance explained", fill = "") +
  scale_fill_viridis_d()
ggsave("figures/area_sensitivity/variance_partition_abundTrend.pdf", units = "in", height = 7, width = 12)  

#### Make a joint model: fixed effects[climate + forest + area_sensitivity] + random slopes/intercepts[species]
# No z score abundance trends, use term for southern limit of wintering grounds

obs_size <- abund_trend %>% 
  group_by(SPEC) %>% 
  summarize(nRoutes = n_distinct(stateroute)) %>%
  arrange(nRoutes) # Smallest = CERW @ 38 routes, next smallest = CAWA @ 48, NOWA @ 54

clim_hab_poptrend_mixedmod <- clim_hab_poptrend_z %>%
  filter(SPEC != "CERW")
#write.csv(clim_hab_poptrend_mixedmod, "\\\\BioArk\\hurlbertlab\\DiCecco\\data\\clim_hab_poptrend_mixedmod.csv", row.names = F)

clim_hab_poptrend_mixedmod <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\clim_hab_poptrend_mixedmod.csv", stringsAsFactors = F)
clim_hab_poptrend_mixedmod <- read.csv("/Volumes/hurlbertlab/DiCecco/data/clim_hab_poptrend_mixedmod.csv", stringsAsFactors = F)



randomslope_add <- lme(abundTrend ~ tmin + tmax + ppt + deltaED + deltaProp + Wintering_Slimit_general, 
                       random = (~deltaProp + deltaED|SPEC), data = clim_hab_poptrend_mixedmod)
randomslope_fullinter <- lme(abundTrend ~ ppt + deltaProp + tmin*deltaED + tmax*deltaED + Wintering_Slimit_general, 
                         random = (~deltaProp + deltaED*tmax|SPEC), data = clim_hab_poptrend_mixedmod)
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
  labs(x = "Trend in Tmin", y = "Abundance trend", color = "Species") 
ggsave("figures/area_sensitivity/mixedmod_randomeffect_tmin.pdf")

# deltaED
ggplot(clim_hab_poptrend_mixedmod, aes(x = deltaED, color = SPEC)) + 
  geom_point(aes(y = sppmean, color = SPEC), alpha = 0.1) +
  geom_line(aes(y = sppmean, color = SPEC), alpha = 0.5, stat="smooth", method = "lm", se = F) + 
  geom_smooth(aes(y = popmean), color = "black", cex = 1, method = "lm", se = F) +
  scale_color_viridis_d() +
  labs(x = "Change in forest edge density", y = "Abundance trend", color = "Species")
ggsave("figures/area_sensitivity/mixedmod_randomeffect_deltaED.pdf")

# Interaction of Tmax and deltaED
clim_hab_poptrend_mixedmod$dEDsign <- ifelse(clim_hab_poptrend_mixedmod$deltaED > 0, "Increase in forest edge density", "Decrease in forest edge density")

theme_set(theme_bw())
ggplot(clim_hab_poptrend_mixedmod, aes(x = tmax, color = SPEC)) + 
  geom_point(aes(y = sppmean, color = SPEC), alpha = 0.1) +
  geom_line(aes(y = sppmean, color = SPEC), alpha = 0.5, stat="smooth", method = "lm", se = F) + 
  geom_smooth(aes(y = popmean), color = "black", cex = 1, method = "lm", se = F) +
  scale_color_viridis_d() +
  labs(x = "Trend in Tmax", y = "Abundance trend", color = "Species") +
  facet_wrap(~dEDsign)
ggsave("figures/area_sensitivity/mixedmod_randomeffect_interaction.pdf", units = "in", height = 6, width = 12)

## Trait analysis
# Predict slope of correlation of spp abundance trends for tmin, tmax, deltaED with three traits
# Mean forest cover of routes observed on
# Mean forest edge density of routes observed on
# Climatic niche breadth

z <- function(x, na.rm = TRUE) {(x - mean(x, na.rm = na.rm))/sd(x, na.rm = na.rm)}

spp_breadths_z <- data.frame(spp_breadths,
                             propFor_z = z(spp_breadths$propFor),
                             ed_z = z(spp_breadths$ed),
                             volume_z = z(spp_breadths$volume))

# Calculate relationship between abundance trend and three env predictors from model for each spp
# Join with spp env traits (z scores)
# Remake these plots!

env_trends_spp <- clim_hab_poptrend_mixedmod %>%
  group_by(aou) %>%
  nest() %>%
  mutate(tmin_mod = map(data, ~{
    df <- .
    mod <- lm(sppmean ~ tmin, data = df)
    tidy(mod)
  }),
  tmax_mod = map(data, ~{
    df <- .
    mod <- lm(sppmean ~ tmax, data = df)
    tidy(mod)
  }),
  ed_mod = map(data, ~{
    df <- .
    mod <- lm(sppmean ~ deltaED, data = df)
    tidy(mod)
  }),
  inter_mod = map(data, ~{
    df <- .
    mod <- lm(sppmean ~ deltaED*tmax, data = df)
  })) %>%
  dplyr::select(-data) %>%
  gather(key = "model", value = "output", 2:4) %>%
  unnest() %>%
  filter(term == "sppmean") %>%
  left_join(spp_breadths_z)

# 3 models, 1 for each slope, slope ~ for_cov + ed + volume

trait_mods <- env_trends_spp %>%
  group_by(model) %>%
  nest() %>%
  mutate(trait_mod = map(data, ~{
    df <- .
    mod <- lm(estimate ~ propFor_z + ed_z + volume_z, data = df)
    tidy(mod)
  })) %>%
  dplyr::select(-data) %>%
  unnest() %>%
  filter(term != "(Intercept)")

ggplot(trait_mods, aes(x = term, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error)) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_grid(~model) +
  labs(x = "", y = "Estimate")
ggsave("figures/area_sensitivity/trait_model_estimates.pdf")

ed_ed <- filter(env_trends_spp, model == "ed_mod") %>% 
  ggplot(aes(x = ed_z, y = estimate)) + 
  geom_smooth(method = "lm", se = F) +
  geom_text(aes(label = SPEC)) +
  labs(x = "Mean edge density", y = "Abundance trend vs. change in ED")

tmin_vol <- filter(env_trends_spp, model == "tmin_mod") %>% 
  ggplot(aes(x = volume_z, y = estimate)) + 
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(x = "Climatic niche breadth", y = "Abundance trend vs. Tmin")

tmin_for <- filter(env_trends_spp, model == "tmin_mod") %>% 
  ggplot(aes(x = propFor_z, y = estimate)) +
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(x = "Mean proportion forest cover", y = "Abundance trend vs. Tmin")

plot_grid(ed_ed, tmin_vol, tmin_for, nrow = 2)
ggsave("figures/area_sensitivity/trait_model_relationships.pdf", units = "in", height = 8, width = 10)

## BRMS model

### Need to use alphadelt bigger than 0.8, see https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
library(brms)
Sys.time()
brm_add <- brm(abundTrend ~ ppt + deltaProp + tmin + deltaED + tmax + Wintering_Slimit_general + (~ppt + deltaProp + tmin + deltaED + tmax + Wintering_Slimit_general|SPEC),
                     data = clim_hab_poptrend_mixedmod,
                     iter = 10000, 
                     cores = 4)
brm_fullinter <- brm(abundTrend ~ ppt + deltaProp + tmin*deltaED + deltaED*tmax + Wintering_Slimit_general + (~ppt + deltaProp + tmin*deltaED + deltaED*tmax + Wintering_Slimit_general|SPEC),
                             data = clim_hab_poptrend_mixedmod,
                     iter = 10000, 
                     cores = 4)
brm_shortinter <- brm(abundTrend ~ tmin*deltaED + deltaED*tmax + Wintering_Slimit_general + (~tmin*deltaED + deltaED*tmax + Wintering_Slimit_general|SPEC),
                     data = clim_hab_poptrend_mixedmod,
                     iter = 10000, 
                     cores = 4)
Sys.time()

params <- summary(brm_add)
fixed <- data.frame(Variable = row.names(params$fixed), params$fixed)

ggplot(filter(fixed, Variable != "Intercept", !(grepl("Wintering", Variable))), aes(x = Variable, y = Estimate)) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = l.95..CI, ymax = u.95..CI)) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  theme_bw()
ggsave("figures/area_sensitivity/brms_fixedeffects_add.pdf", units = "in", height = 6, width = 16)

params <- summary(brm_fullinter)
fixed_fullinter <- data.frame(Variable = row.names(params$fixed), params$fixed)

ggplot(filter(fixed_fullinter, Variable != "Intercept", !(grepl("Wintering", Variable))), aes(x = Variable, y = Estimate)) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = l.95..CI, ymax = u.95..CI)) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  theme_bw()
ggsave("figures/area_sensitivity/brms_fixedeffects_fullinter.pdf", units = "in", height = 6, width = 16)


params <- summary(brm_shortinter)
fixed_shortinter <- data.frame(Variable = row.names(params$fixed), params$fixed)

ggplot(filter(fixed, Variable != "Intercept", !(grepl("Wintering", Variable))), aes(x = Variable, y = Estimate)) +
  geom_point(cex = 2) + geom_errorbar(aes(ymin = l.95..CI, ymax = u.95..CI)) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  theme_bw()
ggsave("figures/area_sensitivity/brms_fixedeffects_shortinter.pdf", units = "in", height = 6, width = 16)

fullinterPlots <- plot(marginal_effects(brm_fullinter, method = "fitted", re_formula = NULL), points = T)
