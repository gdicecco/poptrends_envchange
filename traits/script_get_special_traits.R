## Script: get trait data for specialization in birds (habitat and thermal niche)

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
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

setwd("//BioArk/HurlbertLab/DiCecco/")
correlates <- read.csv("Master_RO_Correlates_20110610.csv", stringsAsFactors = F)
traits <- sppHabit %>%
  left_join(select(correlates, AOU, Tol, OMI), by = c("aou" = "AOU")) %>%
  select(-french_common_name, -spanish_common_name)

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

data(quercus)
head(quercus)
demo('quercus', package='hypervolume')

# Get worldclim at lat/lon occurrences of species

# Test: hairy woodpecker

hairy <- counts.subs %>%
  filter(aou == 3940) %>%
  select(stateroute) %>%
  unique() %>%
  left_join(routes) %>%
  select(longitude, latitude)

climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())

# z-transform climate layers to make axes comparable
climatelayers_ss = climatelayers[[c(1,4,12,15)]]
for (i in 1:nlayers(climatelayers_ss))
   {
  climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
   }
   
climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))

# extract transformed climate values
climate_hairy = raster::extract(climatelayers_ss, hairy)

# unid. error with this function - also happens with iris example and finch example
hairy_hypervol <- hypervolume_gaussian(climate_hairy, name = "hairy", 
                                       kde.bandwidth = estimate_bandwidth(climate_hairy))

# use get_volume() to get volume of niche hypervolume


