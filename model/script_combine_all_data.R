# Build model of habitat fragmentation + climate ~ population trend

######## Reading in and subsetting data ##########
# Population data
## BBS
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

## Abundance indices
setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\data\\")
abundind <- read.csv("BBS_Annual_Indices_Estimates_2015_7-29-2016.csv", stringsAsFactors = F)
poptrends <- read.csv("BBS_Trend_Estimates_2015_7-29-2016.csv", stringsAsFactors = F)
regioncodes <- poptrends %>%
  distinct(Region, Region.Code) %>%
  filter(grepl("[0-9]$", Region.Code)) %>%
  mutate(bcr = as.numeric(substr(Region.Code, 2,3)))

abundind_subs <- abundind %>%
  left_join(regioncodes) %>% 
  mutate(AOU = as.numeric(substr(AOU.Number, 2, 6))) %>%
  filter(bcr %in% bcrs, AOU %in% landbirds$aou, Year >= 1992)

abund_trend <- abundind_subs %>%
  group_by(AOU, bcr) %>%
  nest() %>%
  mutate(lmFit = map(data, ~{
    df <- .
    df.short <- df %>%
      select(Year, Annual.Index) %>%
      unique()
    lm(Annual.Index ~ Year, df.short)
  })) %>%
  mutate(lm_broom = map(lmFit, tidy)) %>%
  mutate(abundTrend = map_dbl(lm_broom, ~{
    df <- .
    slope <- df %>%
      filter(term == "Year") %>%
      select(estimate)
    slope[[1]]
  })) %>%
  mutate(trendPval = map_dbl(lm_broom, ~{
    df <- .
    p <- df %>%
      filter(term == "Year") %>%
      select(p.value)
    p[[1]]
  })) %>%
  filter(abundTrend < 100) # filter out tricolored blackbird

hist(abund_trend$abundTrend, breaks = 100)

# Climate data
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/climate/")
climate_trends <- read.csv("bbs_routes_climate_trends.csv", stringsAsFactors = F)
climate_trends_bcrs <- climate_trends %>%
  left_join(routes) %>%
  group_by(bcr, env) %>%
  summarize(meanTrend = mean(climateTrend))

# Habitat fragmentation data
frags <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\data\\fragmentation_indices_nlcd.csv", stringsAsFactors = F)

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
frags.legend <- bind_rows(frags.92, frags.00s)

frags.filter <- frags.legend %>%
  filter(legend %in% c("Developed, open space", "Developed, low intensity", "Developed, medium intensity", "Developed, high intensity", "Cultivated Crops", "Pasture/hay", "Deciduous forest", "Evergreen forest", "Mixed forest", "Shrub/scrub", "Grassland/herbaceous"))

# Proportion of landscape deltas
frags.01 <- filter(frags.filter, year == 2001)
frags.11 <- filter(frags.filter, year == 2011)
frags.deltas <- frags.11 %>%
  mutate(delta.PL = prop.landscape - frags.01$prop.landscape) %>%
  mutate(delta.ED = edge.density - frags.01$edge.density) %>%
  mutate(delta.PAR = mean.perim.area.ratio - frags.01$mean.perim.area.ratio)

# Trait data
traits <- read.csv("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/traits/spp_traits.csv", stringsAsFactors = F)

############ Build model #########
climate_wide <- climate_trends_bcrs %>%
  spread(key = "env", value = "meanTrend")
  
traits.short <- traits %>%
  select(Common_name, aou, nHabitats1, nHabitats2, volume)

frags.short <- frags.deltas %>%
  select(bcr, class, n.patches, total.area, legend, delta.PL, delta.ED, delta.PAR)

# master data table
clim_hab_poptrend <- abund_trend %>%
  left_join(traits.short, by = c("AOU" = "aou")) %>%
  left_join(frags.short, by = "bcr") %>%
  left_join(climate_wide, by = "bcr")
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/model/")
write.csv(select(clim_hab_poptrend, -data, -lmFit, -lm_broom), "climate_fragmentation_traits_by_species.csv", row.names = F)

mod_clim <- lm(abundTrend ~ ppt + tmax + tmin, clim_hab_poptrend)
plot(clim_hab_poptrend$tmin, clim_hab_poptrend$abundTrend)
abline(mod_clim)
plot(clim_hab_poptrend$ppt, clim_hab_poptrend$abundTrend)
abline(mod_clim)
mod_vol <- lm(abundTrend ~ volume, clim_hab_poptrend)

############ Plots ##############
