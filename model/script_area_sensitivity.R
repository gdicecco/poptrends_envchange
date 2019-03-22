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
            sterr_ed = sd(ED)/sqrt(length(ED)),
            propFor = mean(propForest),
            sterr_for = sd(propForest)/sqrt(length(propForest)),
            patchArea = mean(meanPatchArea),
            sterr_patch = sd(meanPatchArea)/sqrt(length(meanPatchArea)))

ed_breadth <- ggplot(env_breadth, aes(x = reorder(SPEC, ed), y = ed, color = Area_sensitivity)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ed - 1.96*sterr_ed, ymax = ed + 1.96*sterr_ed)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Species", y = "Forest edge density")

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
            sterr_for = sd(propForest)/sqrt(length(propForest)),
            patchArea = mean(meanPatchArea),
            sterr_patch = sd(meanPatchArea)/sqrt(length(meanPatchArea)))

cover_breadth <- ggplot(env_breadth_allroutes, aes(x = reorder(SPEC, propFor), y = propFor, color = Area_sensitivity)) +
  geom_point() + 
  geom_errorbar(aes(ymin = propFor - 1.96*sterr_for, ymax = propFor + 1.96*sterr_for)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Species", y = "Forest cover")

plot_grid(ed_breadth, cover_breadth, nrow = 2)
ggsave("figures/area_sensitivity/forest_breadth.pdf", units = "in", height = 12, width = 15)

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

short_clim <- clim_hab_poptrend_mixedmod %>%
  filter(SPEC %in% obs_size[20:24, ]$SPEC)
randomslope_add <- lme(abundTrend ~ tmin + tmax + ppt + deltaED + deltaProp + Wintering_Slimit_general + Area_sensitivity, 
                       random = (~SPEC|SPEC), data = short_clim)
randomslope_inter <- lme(abundTrend ~ ppt + deltaProp + Wintering_Slimit_general + Area_sensitivity + tmin*deltaED + tmax*deltaED, 
                         random = (~SPEC|SPEC), data = short_clim)

library(brms)
randomslope_add_brm <- brm(abundTrend ~ tmin + tmax + ppt + deltaED + deltaProp + (SPEC|SPEC), 
                       data = clim_hab_poptrend_mixedmod, inits = "0", cores = 4)
randomslope_inter_brm <- brm(abundTrend ~ ppt + deltaProp + Wintering_Slimit_general + Area_sensitivity + tmin*deltaED + tmax*deltaED + (SPEC|SPEC),
                         data = clim_hab_poptrend_mixedmod, inits = "0", cores = 4)

