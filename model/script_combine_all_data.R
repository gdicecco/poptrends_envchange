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
  filter(legend != "Cultivated crops", legend != "Barren land") %>%
  dplyr::select(year, stateroute, edge.density, legend)  %>%
  spread(key = "year", value = "edge.density")
colnames(frags.legend) <-  c("stateroute", "legend", "ED1992", "ED2001", "ED2006", "ED2011")

frag_trends <- frags.legend %>%
  replace_na(list(ED1992 = 0, ED2001 = 0, ED2006 = 0, ED2011 = 0)) %>%
  mutate(deltaED = ED2011 - ED1992)

# Trait data
traits <- read.csv("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/traits/spp_traits.csv", stringsAsFactors = F)

############ Build model #########
  
traits.short <- traits %>%
  dplyr::select(Common_name, aou, nHabitats1, nHabitats2, volume)

# master data table
clim_hab_poptrend <- abund_trend %>%
  left_join(traits.short) %>%
  left_join(frag_trends, by = "stateroute") %>%
  left_join(climate_wide, by = "stateroute")
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/model/")
write.csv(clim_hab_poptrend, "climate_fragmentation_traits_by_species.csv", row.names = F)

## Start with 10 most abundant species
abund_spp <- counts.subs %>% 
  group_by(aou) %>%
  summarize(spptotal = sum(speciestotal)) %>%
  arrange(desc(spptotal)) %>%
  slice(1:20) %>%
  left_join(species)

############ Plots ##############
