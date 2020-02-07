### Calculate forest specialization of species using mean forest cover of BBS routes
### Create traits table used for trait models in model/script_main_analyses.R

library(tidyverse)

## Make master data table

abund_trend <- read.csv("model/BBS_abundance_trends.csv", stringsAsFactors = F)

env_change <- read.csv("model/bbs_route_env_change.csv", stringsAsFactors = F)

fourletter_codes <- read.csv("traits/four_letter_codes_aou.csv", stringsAsFactors = F)

clim_hab_poptrend <- abund_trend %>%
  left_join(env_change, by = "stateroute") %>%
  filter(!is.na(propForest))

## Read in forest cover data for BBS routes

frags <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_nlcd_simplified.csv", stringsAsFactors = F)
frags_ca <- read.csv("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_canada.csv", stringsAsFactors = F) %>%
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

## Read in species trait data from Hurlbert & White 2007

correlates <- read.csv("//BioArk/HurlbertLab/DiCecco/Data/Master_RO_Correlates_20110610.csv", stringsAsFactors = F)

#### Species traits ####
# Breadth of forest ED and forest cover for area and non-area sensitive species

env_breadth <- clim_hab_poptrend %>%
  left_join(fourletter_codes) %>%
  dplyr::group_by(SPEC, aou) %>%
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
  left_join(forest_allroutes, by = "stateroute") %>%
  left_join(dplyr::select(env_change, stateroute, year, tmax, tmin), by = "stateroute") %>%
  filter(!is.na(propForest))

env_breadth_allroutes <- clim_hab_pop_allroutes %>%
  left_join(fourletter_codes) %>%
  dplyr::group_by(SPEC, aou) %>%
  summarize(propFor = mean(propForest),
            min_for = quantile(propForest, c(0.05))[[1]],
            max_for = quantile(propForest, c(0.95))[[1]],
            patchArea = mean(meanPatchArea),
            std_patch = sd(meanPatchArea))

spp_breadths <- env_breadth_allroutes %>%
  dplyr::select(SPEC, aou, propFor, min_for, max_for) %>%
  left_join(dplyr::select(env_breadth, SPEC, aou, ed, std_ed))

# Species table for MS supplement

temp_range <- read.csv("traits/bbs_aou_temp_range.csv", stringsAsFactors = F)

spp_table_traits <- correlates %>%
  dplyr::select(AOU, CommonName, migclass, Foraging, Brange_Area_km2) %>%
  left_join(spp_breadths, by = c("AOU" = "aou")) %>%
  left_join(temp_range, by = c("AOU" = "aou")) %>%
  filter(!is.na(SPEC), SPEC != "CERW") %>%
  ungroup() %>%
  dplyr::select(AOU, CommonName, SPEC, migclass, Foraging, propFor, temp_range, Brange_Area_km2)

# write.csv(spp_table_traits, "traits/forest_spp_traits_MS.csv", row.names = F)

