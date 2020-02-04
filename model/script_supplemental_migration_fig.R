# Supplemental figure: migration connectivity

library(tidyverse)
library(sf)
library(tmap)

abund_trend <- read.csv("model/BBS_abundance_trends.csv", stringsAsFactors = F)

#### Migration connectivity figure ####

# For Rose-breasted grosbeak, Black-throated blue warbler, Wood thrush, Cedar waxwing
# Make maps of abundance trends (smoothed?), showing wintering grounds (connect with lines in PPT?)

# Read in range maps

americas_map <- world %>%
  filter(region_un == "Americas")

seasonal_code <- data.frame(SEASONAL = 1:5,
                            Legend = c("Resident", "Breeding", "Non-breeding", "Passage", "Seasonal occurrence uncertain"),
                            stringsAsFactors = F)

rogr_range <- read_sf("\\\\BioArk//HurlbertLab//GIS//birds//All//All//Pheucticus_ludovicianus_22723813.shp") %>%
  filter(PRESENCE < 5, SEASONAL < 4) %>%
  left_join(seasonal_code)

rogr_abund <- abund_trend %>%
  left_join(routes) %>%
  filter(Common_name == "Rose-breasted grosbeak") %>%
  st_as_sf(coords = c("longitude", "latitude"))

rogr_map <- tm_shape(rogr_range) + tm_polygons(col = "Legend", palette = "-Greys", title = " ") + 
  tm_layout(legend.text.size = 2) +
  tm_shape(americas_map) + tm_borders() +
  tm_shape(rogr_abund) + tm_dots(col = "abundTrend", size = 0.1, palette = "-RdBu", title = "Abundance trend") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.25, main.title = "A. Rose-breasted grosbeak")

btbl_range <- read_sf("\\\\BioArk//HurlbertLab//GIS//birds//All//All//Dendroica_caerulescens_22721673.shp") %>%
  filter(PRESENCE < 5, SEASONAL < 4) %>%
  left_join(seasonal_code)

btbl_abund <- abund_trend %>%
  left_join(routes) %>%
  filter(Common_name == "Black-throated blue warbler") %>%
  st_as_sf(coords = c("longitude", "latitude"))

btbl_map <- tm_shape(btbl_range) + tm_polygons(col = "Legend", palette = "-Greys", title = " ", legend.show = F) + 
  tm_layout(legend.text.size = 2) +
  tm_shape(americas_map) + tm_borders() +
  tm_shape(btbl_abund) + tm_dots(col = "abundTrend", size = 0.1, palette = "-RdBu", title = "Abundance trend") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.25, main.title = "B. Black-throated blue warbler",
            legend.position = c(0.65, 0.4))

woth_range <- read_sf("\\\\BioArk//HurlbertLab//GIS//birds//All//All//Hylocichla_mustelina_22708670.shp") %>%
  filter(PRESENCE < 5, SEASONAL < 4) %>%
  left_join(seasonal_code)

woth_abund <- abund_trend %>%
  left_join(routes) %>%
  filter(Common_name == "Wood thrush") %>%
  st_as_sf(coords = c("longitude", "latitude"))

woth_map <- tm_shape(woth_range) + tm_polygons(col = "Legend", palette = "-Greys", title = " ", legend.show = F) + 
  tm_layout(legend.text.size = 2) +
  tm_shape(americas_map) + tm_borders() +
  tm_shape(woth_abund) + tm_dots(col = "abundTrend", size = 0.1, palette = "-RdBu", title = "Abundance trend") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.25, main.title = "C. Wood thrush",
            legend.position = c(0.68, 0.3))

cewa_range <- read_sf("\\\\BioArk//HurlbertLab//GIS//birds//All//All//Bombycilla_cedrorum_22708153.shp") %>%
  filter(PRESENCE < 5, SEASONAL < 4) %>%
  left_join(seasonal_code)

cewa_range$Legend_relevel <-  fct_relevel(cewa_range$Legend, c("Breeding", "Resident", "Non-breeding"))

cewa_abund <- abund_trend %>%
  left_join(routes) %>%
  filter(Common_name == "Cedar waxwing") %>%
  st_as_sf(coords = c("longitude", "latitude"))

cewa_map <- tm_shape(cewa_range) + tm_polygons(col = "Legend_relevel", palette = "-Greys", title = " ") + 
  tm_layout(legend.text.size = 2) +
  tm_shape(americas_map) + tm_borders() +
  tm_shape(cewa_abund) + tm_dots(col = "abundTrend", size = 0.08, palette = "-RdBu", title = "Abundance trend") +
  tm_layout(legend.text.size = 0.8, legend.title.size = 1.25, main.title = "D. Cedar waxing",
            legend.position = c("left", "bottom"))

mig_con_fig <- tmap_arrange(rogr_map, btbl_map, woth_map, cewa_map, ncol = 2)
tmap_save(mig_con_fig, "figures/main_analysis_figs/migration_connectivity_fig.pdf", height = 9, width = 9, unit = "in")
