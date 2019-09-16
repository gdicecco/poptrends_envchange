# Explore abundance trends

library(tidyverse)
library(sp)
library(tmap)
library(sf)

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
setwd("/Volumes/hurlbertlab/DiCecco/data/")
abund_trend <- read.csv("BBS_abundance_trends.csv", stringsAsFactors = F) %>%
  filter(aou %in% area_aous$aou) %>% ## Only ones lost from area_aous are nocturnal: owls, nightjars (Caprimulgus)
  left_join(area_aous) %>%
  filter(sporder != "Accipitriformes" & sporder != "Strigiformes" & sporder != "Caprimulgiformes") ## Remove hawks

length(abund_trend$trendPval[abund_trend$trendPval < 0.05])/nrow(abund_trend) # 44%
length(abund_trend$trendPval[abund_trend$trendPval < 0.01])/nrow(abund_trend) # 35%
length(abund_trend$trendPval[abund_trend$trendPval < 0.0001])/nrow(abund_trend) # 21%

abund_trend_sig <- abund_trend %>%
  mutate(abund_class = case_when(trendPval >= 0.05 ~ "p > 0.05",
                                 trendPval >= 0.01 & trendPval < 0.05 ~ "0.01 < p < 0.05",
                                 trendPval >= 0.001 & trendPval < 0.01 ~ "0.001 < p < 0.01",
                                 trendPval < 0.001 ~ "p < 0.001"),
         obs_class = case_when(obsPval >= 0.05 ~ "p > 0.05",
                               obsPval >= 0.01 & obsPval < 0.05 ~ "0.01 < p < 0.05",
                               obsPval >= 0.001 & obsPval < 0.01 ~ "0.001 < p < 0.01",
                               obsPval < 0.001 ~ "p < 0.001"))

setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/")

# Histogram of abundance trends colored by p value significance level

ggplot(abund_trend_sig, aes(x = abundTrend, 
                            fill = fct_relevel(abund_class, c("p < 0.001",
                                                              "0.001 < p < 0.01", 
                                                              "0.01 < p < 0.05", 
                                                              "p > 0.05")))) + 
  geom_histogram(bins = 20) + 
  labs(fill = "", x = "Abundance trend", y = "Count") + 
  scale_fill_viridis_d()
ggsave("figures/area_sensitivity/hist_abundTrend_sig.pdf")

# Individual maps for species of abundance trends in range

us.proj <- readOGR("\\\\BioArk/hurlbertlab/DiCecco/nlcd_frag_proj_shapefiles/BCRs_contiguous_us/BCRs_contiguous_us.shp")
states <- us.proj[, -(1:2)]

abund_trend_map <- abund_trend %>%
  left_join(routes) %>%
  st_as_sf(coords = c("longitude", "latitude"))

pdf("figures/area_sensitivity/species_abundance_maps.pdf")
for(spp in unique(abund_trend$aou)) {
  df <- filter(abund_trend_map, aou == spp)
  spec <- fourletter_codes$SPEC[fourletter_codes$aou == spp]
  
  print(tm_shape(states) + tm_polygons(col = 'white') + tm_shape(df) + tm_dots(col = "abundTrend", palette = "RdBu", size = 1) +
          tm_layout(title = spec))
}
dev.off()

# Individual species maps - average by state

states_sf <- read_sf("\\\\BioArk\\Hurlbertlab\\GIS\\geography\\ne_50m_admin_1_states_provinces_lakes\\ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA") %>%
  dplyr::select(adm1_code) %>%
  filter(adm1_code != "USA-3563" & adm1_code != "USA-3517")
crs <- st_crs(states_sf)

pdf("figures/area_sensitivity/species_abundance_maps_states.pdf")
for(spp in unique(abund_trend$aou)) {
  
  spec <- fourletter_codes$SPEC[fourletter_codes$aou == spp]
  
  abund_states <- abund_trend_map %>%
    filter(aou == spp) %>%
    st_set_crs(crs) %>%
    st_join(states_sf) %>%
    group_by(adm1_code) %>%
    summarize(meanTrend = mean(abundTrend, na.rm = T))
  
  abund_states_df <- data.frame(adm1_code = abund_states$adm1_code, meanTrend = abund_states$meanTrend)
  
  abund_states_polygon <- states_sf %>%
    left_join(abund_states_df)
  
  print(tm_shape(abund_states_polygon) + tm_borders() + 
          tm_fill(col= "meanTrend", palette = "RdBu") + tm_layout(title = spec))
  
}
dev.off()

# For one species with not that many routes - show counts vs. time

bhvi <- filter(abund_trend, aou == 6290, nObs > 9)

pdf("figures/area_sensitivity/abundance_trends_bhvi.pdf")
for(strte in bhvi$stateroute) {
  df <- filter(bhvi, stateroute == strte)
  glm <- df$lmFit[[1]]
  df.count <- df$data[[1]]
  df.count$predict <- predict(glm, newdata = df.count, type = "response")
  plot <- ggplot(df.count, aes(x = year, y = speciestotal)) + geom_point() + 
    geom_smooth(method = "glm", se = F, method.args = list(family = "poisson"))
  
  print(plot)
  
}
dev.off()

# Histogram of first year observer effect colored by p value significance level

ggplot(abund_trend_sig, aes(x = obsTrend, 
                            fill = fct_relevel(obs_class, c("p < 0.001",
                                                            "0.001 < p < 0.01", 
                                                            "0.01 < p < 0.05", 
                                                            "p > 0.05")))) + 
  geom_histogram(bins = 20) + 
  labs(fill = "", x = "First year observer effect", y = "Count") + 
  scale_fill_viridis_d()
ggsave("figures/area_sensitivity/hist_obsEffect_sig.pdf")

# Which species most commonly show observer effects?

presences <- counts.subs %>%
  group_by(aou) %>%
  summarize(range_size = n_distinct(stateroute))

obs_effects <- abund_trend_sig %>% 
  filter(!is.na(obs_class), obsPval < 0.05) %>%
  group_by(aou) %>%
  count() %>%
  left_join(presences) %>%
  left_join(fourletter_codes) %>%
  mutate(obs_prop = n/range_size) %>%
  arrange(desc(obs_prop))
write.csv(obs_effects, "model/obs_effect_rankings.csv", row.names = F)
