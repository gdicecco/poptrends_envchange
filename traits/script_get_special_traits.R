## Script: get trait data for specialization in birds (habitat and thermal niche)

library(dplyr)
library(stringr)
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
  filter(aou %in% unique(counts.subs$aou))

# Which species didn't match up
checklist.unid <- checklist.subs[is.na(checklist.subs$SISRecID), ] #NAs for SISRecID
checklist.nas <- checklist.unid[!grepl("unid.", checklist.unid$english_common_name), ] # Omit the ones that are NA because they are unid.

# Manually entered missing woodpeckers, warblers, sparrows
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/traits/")
missingspp <- read.csv("spp_missing_ids.csv", stringsAsFactors = F)

## Index from 0-1 for habitats - logit transform

# Get habitat data
IUCNids <- na.omit(unique(checklist.subs$SISRecID))
finescale_habitats <- matrix(nrow = length(IUCNids), ncol = 2) 

for(i in 1:length(IUCNids)) {
  id <- IUCNids[i]
  finescale_habitats[i,1] <- id
  habitat <- birdlife_habitat(id)
  finescale_habitats[i,2] <- length(habitat$id)
}
habitats <- data.frame(SISRecID = finescale_habitats[,1], nHabitats = finescale_habitats[,2]) #38 IDs are zero

habitats$index <- 1 - (habitats$nHabitats - 1)/(max(habitats$nHabitats) - 1)
hist(habitats$index)

sppHabit <- checklist.subs %>%
  left_join(habitats) %>%
  arrange(nHabitats)

## Thermal niche 
# Start with Tol - high tolerance = broad niche

setwd("//BioArk/HurlbertLab/DiCecco/")
correlates <- read.csv("Master_RO_Correlates_20110610.csv", stringsAsFactors = F)
traits <- sppHabit %>%
  left_join(select(correlates, AOU, Tol, OMI), by = c("aou" = "AOU")) %>%
  select(-french_common_name, -spanish_common_name)

## Compare thermal niche breadth and habitat specialization
library(ggplot2)
theme_set(theme_bw())

traitplot <- ggplot(traits, aes(x = index, y = Tol)) + geom_point() 
traitplot + xlab("Rel. habitat specialization") + ylab("Tolerance") + geom_smooth(method = "lm", col = "blue", se = F)
# Increased habitat specialization is correlated (weakly) with decreased thermal niche breadth


# Percentile - specialization

