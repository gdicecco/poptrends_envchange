## Script: get trait data for specialization in birds (habitat and thermal niche)

library(dplyr)
library(stringr)

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
  filter(sporder == "Piciformes" | sporder == "Passeriformes" | sporder == "Psittaciformes" | sporder == "Coraciiformes" |
         sporder == "Trogoniformes" | sporder == "Apodiformes"| sporder == "Columbiformes" | sporder == "Cuculiformes" |
         sporder == "Galliformes")

# Observations with RT=1, in BCRs of interest, 1990-present, diurnal land birds
counts.subs <- counts %>%
  filter(stateroute %in% routes.short$stateroute) %>%
  filter(year >= 1990) %>%
  filter(aou %in% landbirds$aou)

# Number of species: 342
length(unique(counts.subs$aou))

## Get number of habitats
setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\LTER_birdabund_seasonal\\")

##BirdLife checklist (ID numbers)
checklist <- read.csv("BirdLife_Checklist_V_9.1.csv", header = TRUE, stringsAsFactors = F)
checklist.subs <- checklist[,c(4,5,9,10,13)] #just SISRecID, scientific name, synonyms, alt common names, common name
  
checklist.subs <- checklist %>%
  select(Common_name, Scientific_name, Synonyms, Alt_common_names, SISRecID) %>%
  mutate(genus = word(Scientific_name, 1),
         species = word(Scientific_name, 2),
         alt_genus = word(Synonyms, 1),
         alt_spp = word(Synonyms, 2)) %>%
  right_join(landbirds, by = c("genus", "species")) %>% 
  filter(aou %in% unique(counts.subs$aou))

checklist.unid <- checklist.subs[is.na(checklist.subs$SISRecID), ]
checklist.nas <- checklist.unid[!grepl("unid.", checklist.unid$english_common_name), ]

# Manually entered woodpeckers, missing warblers, missing sparrows
missingspp <- read.csv("spp_missing_ids.csv", stringsAsFactors = F)

## Index from 0-1 for habitats - logit transform

## Thermal niche 

## Thermal niche index 0-1

