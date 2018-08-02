## Script: get trait data for specialization in birds (habitat and thermal niche)

library(tidyverse)
library(traits)

## List of species observed in BCRs of interest during time window (1990-present)
routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")

bcrs <- c(9, 12, 13, 14, 18, 19, 23, 27, 29) # BCRs of interest

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, bcr %in% bcrs, select = c("statenum","stateroute", "year", "latitude", "longitude", "bcr"))
counts$stateroute <- counts$statenum*1000 + counts$route

landbirds <- species %>%
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010)

# Observations with RT=1, in BCRs of interest, 1990-present, diurnal land birds
counts.subs <- counts %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  filter(year >= 1990) %>%
  filter(aou %in% landbirds$aou)

# Number of species: 367
length(unique(counts.subs$aou))

## Get number of habitats
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\LTER_birdabund_seasonal\\")

## BirdLife checklist (ID numbers)
checklist <- read.csv("BirdLife_Checklist_V_9.1.csv", header = TRUE, stringsAsFactors = F)
  
checklist.subs <- checklist %>%
  select(Common_name, Scientific_name, Synonyms, Alt_common_names, SISRecID) %>%
  mutate(genus = word(Scientific_name, 1),
         species = word(Scientific_name, 2)) %>%
  right_join(landbirds, by = c("genus", "species")) %>% 
  replace_na(list(Common_name = "unknown")) %>%
  filter(aou %in% unique(counts.subs$aou), Common_name != "")

# Check which species didn't match up
checklist.unid <- checklist.subs[is.na(checklist.subs$SISRecID), ] #NAs for SISRecID
checklist.nas <- checklist.unid[!grepl("unid.", checklist.unid$english_common_name), ] # Omit the ones that are NA because they are unid.

# Manually entered species (taxonomic changes in spp_taxon_changes.csv)
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/traits/")
missingspp <- read.csv("spp_missing_ids.csv", stringsAsFactors = F)

checklist.subs$SISRecID[match(missingspp$aou, checklist.subs$aou)] <- missingspp$SISRecID

# Get habitat data
IUCNids <- na.omit(unique(checklist.subs$SISRecID))
finescale_habitats <- matrix(nrow = 1, ncol = 5) 
colnames(finescale_habitats) <- c("id", "habitat1", "habitat2", "importance", "occurrence")

for(i in 1:length(IUCNids)) {
  id <- IUCNids[i]
  habitat <- birdlife_habitat(id)
  colnames(habitat) <- c("id", "habitat1", "habitat2", "importance", "occurrence")
  finescale_habitats <- rbind(finescale_habitats, habitat)
}

habitats <- finescale_habitats %>%
  filter(occurrence == "breeding" | occurrence == "resident") %>%
  group_by(id) %>%
  summarize(nHabitats1 = length(unique(habitat1)), nHabitats2 = n())

importance <- finescale_habitats %>%
  filter(occurrence == "breeding" | occurrence == "resident") %>%
  group_by(id, importance) %>%
  summarize(nHabitats = n()) %>%
  arrange(id) %>%
  spread(importance, nHabitats) %>%
  replace_na(list(major = 0, marginal = 0, suitable = 0)) %>%
  mutate(weighted = (marginal*1 + suitable*2 + major*3)/sum(marginal, suitable, major)) # This is super arbitrary

# Compare major vs suitable habitats
theme_set(theme_bw())
implot <- ggplot(importance, aes(x = major, y = suitable)) + geom_point() + geom_smooth(method = "lm", se = F)
implot + ylab("Suitable habitats") + xlab("Major habitats")

# Compare level 1 vs level 2 habitats
habplot <- ggplot(habitats, aes(x = nHabitats1, y = nHabitats2)) + geom_point() + geom_smooth(method = "lm", se = F)
habplot + xlab("Level 1 Habitats") + ylab("Level 2 Habitats")

habitats$index <- 1 - (habitats$nHabitats2 - 1)/(max(habitats$nHabitats2) - 1)
hist(habitats$index)

sppHabit <- checklist.subs %>%
  left_join(habitats, by = c("SISRecID" = "id")) %>%
  left_join(importance, by = c("SISRecID" = "id"))

sppHabit.missing <- sppHabit %>%
  filter(is.na(nHabitats1))

# Compare level 1 to level 2 habitats
unique(finescale_habitats$habitat1)
unique(finescale_habitats$habitat2)
# Level 1 is similar to 92-01 classifications
# Level 2 is more comparable to 01-11 classifications, but no distinction in urban areas
# Forest distinctions from IUCN are boreal/temperate/subtropical/tropical vs. coniferous and deciduous in 01-11 NLCD

## Thermal niche 
# Start with Tol - high tolerance = broad niche

setwd("//BioArk/HurlbertLab/DiCecco/Data/")
correlates <- read.csv("Master_RO_Correlates_20110610.csv", stringsAsFactors = F)
traits <- sppHabit %>%
  left_join(dplyr::select(correlates, AOU, Tol, OMI), by = c("aou" = "AOU")) %>%
  dplyr::select(-french_common_name, -spanish_common_name)

## Compare thermal niche breadth and habitat specialization

traitplot <- ggplot(traits, aes(x = index, y = Tol)) + geom_point() 
traitplot + xlab("Rel. habitat specialization") + ylab("Tolerance") + geom_smooth(method = "lm", col = "blue", se = F)

traitplot.raw <- ggplot(traits, aes(x = weighted, y = Tol)) + geom_point() 
traitplot.raw + xlab("No. habitats used") + ylab("Tolerance") + geom_smooth(method = "lm", col = "blue", se = F)
# Increased habitat specialization is correlated (weakly) with decreased thermal niche breadth
# Changes when weighted average is used 

# Hypervolume to quantify thermal niche width

library(hypervolume)
require(raster)
require(maps)

# WorldClim
climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
climatelayers_ss = climatelayers[[c(10, 18)]]

# z-transform climate layers to make axes comparable
for (i in 1:nlayers(climatelayers_ss))
{
  climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
}

climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))

# Elevation
elev <- raster("\\\\BioArk\\HurlbertLab\\GIS\\DEM\\USA1_msk_alt.grd")
elev_z <- (elev - cellStats(elev, 'mean'))/cellStats(elev, 'sd')

# NDVI data
### CSV obtained using script  https://github.com/ethanwhite/bbs-forecasting/blob/master/R/get_ndvi_data.R 
gimms_ndvi = read.csv("https://raw.githubusercontent.com/hurlbertlab/Biotic-Interactions/master/ENV%20DATA/gimms_ndvi_bbs_data.csv?token=AdzxAN_oTo30YIQcPTsjtzuK0GbPx0Lfks5baevMwA%3D%3D", header = TRUE)
gimms_agg = gimms_ndvi %>% filter(month == c("may", "jun", "jul")) %>% 
  group_by(site_id)  %>%  summarise(ndvi.mean=mean(ndvi))
gimms_agg$stateroute = gimms_agg$site_id
ndvi = gimms_agg[,c("stateroute", "ndvi.mean")] %>%
  mutate(ndvi.z = (ndvi.mean - mean(ndvi.mean))/sd(ndvi.mean)) %>%
  left_join(routes) %>%
  dplyr::select(stateroute, longitude, latitude, ndvi.z)

## Species observed on 10 or more routes
birds.subs <- counts.subs %>%
  group_by(aou) %>%
  summarize(nRoutes = length(unique(stateroute))) %>%
  filter(nRoutes >= 10) %>%
  left_join(counts.subs)

## Takes ~ 2 hours
bird.coords <- birds.subs %>%
  group_by(aou) %>%
  dplyr::select(stateroute) %>%
  unique() %>%
  left_join(routes) %>%
  dplyr::select(stateroute, longitude, latitude) %>%
  nest() %>%
  mutate(volume = purrr::map_dbl(data, ~{
    df <- .
    climate <- raster::extract(climatelayers_ss_cropped, dplyr::select(df, longitude, latitude))
    elev <- raster::extract(elev_z, dplyr::select(df, longitude, latitude))
    env <- ndvi %>%
      right_join(df) %>%
      dplyr::select(ndvi.z) %>%
      cbind(climate, elev)
    hypervol <- hypervolume_gaussian(climate)
    get_volume(hypervol)
  }))

birds.vol <- dplyr::select(bird.coords, aou, volume)

traits <- traits %>% left_join(birds.vol)
write.csv(traits, "C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/traits/spp_traits.csv", row.names = F)

# Compare hypervolume to tolerance
ggplot(traits, aes(x = Tol, y = volume)) + geom_point() + xlab("Tolerance") + ylab("Hypervolume") + geom_smooth(method = "lm", se = F)
