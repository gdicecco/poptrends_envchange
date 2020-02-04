### Fit species-specific models using area-sensitivity classification

#### Libraries ####

library(tidyverse)
library(raster)
library(rgdal)
library(purrr)
library(broom)
library(tmap)
library(cowplot)
library(sf)
library(forcats)
library(grid)
library(spdep)

##### Analysis #####

## Make master data table

abund_trend <- read.csv("model/BBS_abundance_trends.csv", stringsAsFactors = F)

env_change <- read.csv("model/bbs_route_env_change.csv", stringsAsFactors = F)

fourletter_codes <- read.csv("traits/four_letter_codes_aou.csv", stringsAsFactors = F)

spp_table_traits <- read.csv("traits/forest_spp_traits_MS.csv", stringsAsFactors = F)

clim_hab_poptrend <- abund_trend %>%
  left_join(env_change, by = "stateroute") %>%
  filter(!is.na(propForest))

#### Individual models ####
# of climate + frag + loss + climate:frag + climate:loss for species

z <- function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T)}

forest_z <- data.frame(stateroute = env_change$stateroute,
                       ED = z(env_change$ED),
                       propForest = z(env_change$propForest),
                       deltaED = z(env_change$deltaED), 
                       deltaProp = z(env_change$deltaProp))

climate_wide_z <- env_change %>%
  dplyr::select(stateroute, tmin, tmax) %>%
  mutate_at(.vars = c("tmin", "tmax"),
            .funs = z)

clim_hab_poptrend_z <- abund_trend %>%
  left_join(forest_z, by = "stateroute") %>%
  left_join(climate_wide_z, by = "stateroute") %>%
  filter(!is.na(propForest)) %>%
  left_join(fourletter_codes)

# Spatial CAR models by species

model_fits <- clim_hab_poptrend_z %>%
  group_by(aou, SPEC) %>%
  nest() %>%
  mutate(lmFit = purrr::map(data, ~{
    df <- .

    # get species weights matrix

    # Coordinates of BBS routes

    coords_bbs <- df %>%
      ungroup() %>%
      dplyr::select(stateroute) %>%
      distinct() %>%
      left_join(routes) %>%
      dplyr::select(stateroute, longitude, latitude)

    coords_mat <- as.matrix(coords_bbs[, -1])

    # Calculate nearest neighbors and make spatial neighborhood
    k0 <- knearneigh(coords_mat, longlat=T, k=1)
    k1 <- knn2nb(k0)

    # Find maximum neighbor distance and use this distance to define neighborhoods
    max.d <- max(unlist(nbdists(k1, coords_mat, longlat=T)))
    nb0.birds <- dnearneigh(coords_mat, 0, max.d, longlat=T)
    plot(nb0.birds, coords_mat)

    # Using a distance threshold of 100km
    nb1.birds <- dnearneigh(coords_mat,1,100,longlat=T)
    plot(nb1.birds, coords_mat)

    # Create spatial weights based on linear distance decay
    glist.birds <- nbdists(nb0.birds, coords_mat, longlat=T)
    glist.birds <- lapply(glist.birds, function(x) 1-(x/ceiling(max.d))) # Round up max.d so that all point get weights
    wt.birds <- nb2listw(nb0.birds, style='B', glist=glist.birds)

    spautolm(abundTrend ~ tmax*deltaED + tmin*deltaED + deltaProp, data = df, wt.birds, na.action = na.fail, family = "CAR")
  })) %>%
  mutate(pred = purrr::map(lmFit, ~{
    mod <- .
    data.frame(pred.vals = mod$fit$fitted.values)
  })) %>%
  mutate(comb = map2(data, pred, ~bind_cols(.x, .y))) %>%
  mutate(r2 = map_dbl(comb, ~{
    df <- .
    summary(lm(abundTrend ~ pred.vals, df))$r.squared
  })) %>%
  mutate(nObs = map_dbl(data, ~{
    df <- .
    nrow(df)
  }))  %>%
  mutate(lm_broom = purrr::map(lmFit,  ~{
    mod <- .
    sum <- summary(mod)
    df <- as.data.frame(sum$Coef)
    df$term = row.names(df)
    df
  })) %>%
  dplyr::select(aou, SPEC, nObs, lm_broom, r2) %>%
  unnest() %>%
  filter(nObs > 40) %>% # 1 spp at 38, only one below 40
  filter(term != "(Intercept)") %>%
  rename(std.error = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`) %>%
  mutate(sig = case_when(p.value < 0.05 ~ "p < 0.05",
                         TRUE ~ "p > 0.05"))

range(model_fits$nObs)
# 67 species

# write.csv(model_fits, "model/individual_species_model_tables.csv", row.names = F)

model_fits <- read.csv("model/individual_species_model_tables.csv", stringsAsFactors = F)

# Figure: distributions of effect sizes by species

model_fits <- model_fits %>%
  mutate(sig = case_when(p.value < 0.05 & p.value >= 0.01 ~ "0.01 < p < 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "0.001 < p < 0.01",
                         p.value < 0.001 ~ "p < 0.001",
                         TRUE ~ "p > 0.05")) %>%
  mutate(sig = fct_relevel(sig, c("p < 0.001", "0.001 < p < 0.01", "0.01 < p < 0.05", "p > 0.05"))) %>%
  mutate(dir = ifelse(Estimate < 0, "Negative effect", "Positive effect"))

density_plot <- function(df, variable, label) {
  ggplot(filter(df, term == variable), aes(x = Estimate, fill = sig)) + 
    geom_histogram(bins = 20) + 
    geom_vline(aes(xintercept = mean(Estimate)), lty = 2, cex = 1) +
    labs(x = label, y = "Species", fill = "") +
    scale_fill_manual(values = c("p < 0.001" = "#2166AC", 
                                 "0.001 < p < 0.01" = "#67A9CF",
                                 "0.01 < p < 0.05" = "#D1E5F0",
                                 "p > 0.05" = "gray"))
}

deltaED <- density_plot(model_fits, "deltaED", "Change in edge density")
deltaProp <- density_plot(model_fits, "deltaProp", "Change in forest cover")
tmin <- density_plot(model_fits, "tmin", "Trend in Tmin")
tmax <- density_plot(model_fits, "tmax", "Trend in Tmax")
tmaxED <- density_plot(model_fits, "tmax:deltaED", "Tmax:change in ED")
tminED <- density_plot(model_fits, "deltaED:tmin", "Tmin:change in ED")

legend <- get_legend(deltaED)

grid_effects <- plot_grid(deltaED + theme(legend.position = "none"), 
                          deltaProp + theme(legend.position = "none") + ylab(" "), 
          tmin + theme(legend.position = "none"), 
          tmax + theme(legend.position = "none") + ylab(" "), 
          tmaxED + theme(legend.position = "none") + scale_x_continuous(breaks = c(-0.01, 0, 0.01)), 
          tminED + theme(legend.position = "none") + ylab(" "),
          nrow = 3,
          labels = c("A", "B", "C", "D", "E", "F"))
plot_grid(grid_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/model_effects_distributions.pdf", units = "in", height = 9, width = 10)

## Supplemental figure: p value ranges for each variable

pval_plot <- function(df, variable, label) {
  ggplot(filter(df, term == variable), aes(x = p.value, fill = dir)) + 
    geom_histogram(bins = 20) + 
    geom_vline(aes(xintercept = mean(p.value)), lty = 2, cex = 1) + 
    labs(x = label, y = "Count", fill = "")  +
    scale_x_log10() + scale_fill_manual(values = c("dodgerblue3", "firebrick2"))

}

options(scipen = 999)
deltaED <- pval_plot(model_fits, "deltaED", "Change in edge density") + scale_x_log10(breaks = c(0.00001, 0.001, 1),
                                                                labels = c(0.00001, 0.001, 1))
deltaProp <- pval_plot(model_fits, "deltaProp", "Change in forest cover")
tmin <- pval_plot(model_fits, "tmin", "Trend in Tmin")
tmax <- pval_plot(model_fits, "tmax", "Trend in Tmax")
tmaxED <- pval_plot(model_fits, "tmax:deltaED", "Tmax:change in ED")
tminED <- pval_plot(model_fits, "deltaED:tmin", "Tmin:change in ED")

legend <- get_legend(tmin)

grid_effects <- plot_grid(deltaED + theme(legend.position = "none"), 
                          deltaProp + ylab(" ") + theme(legend.position = "none"), 
                          tmin + theme(legend.position = "none"), 
                          tmax + ylab(" ") + theme(legend.position = "none"), 
                          tmaxED + theme(legend.position = "none"), 
                          tminED + ylab(" ") + theme(legend.position = "none"),
                          nrow = 3,
                          labels = c("A", "B", "C", "D", "E", "F"))
plot_grid(grid_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/model_pvals_distributions.pdf", units = "in", height = 9, width = 10)

#### Trait models ####

spp_traits <- spp_table_traits %>%
  left_join(model_fits) %>%
  mutate(foraging_alph = case_when(Foraging == "foliage glean" ~ paste("a", Foraging),
                                   TRUE ~ Foraging)) %>%
  group_by(term) %>%
  nest() %>%
  filter(!is.na(term)) %>%
  mutate(trait_mod = purrr::map(data, ~{
    df <- .
    lm(Estimate ~  temp_range + propFor + log10(Brange_Area_km2) + migclass + foraging_alph, data = df)
  }),
  tidy = purrr::map(trait_mod, ~{
    mod <- .
    tidy(mod)
  }),
  r2 = map_dbl(trait_mod, ~{
    mod <- .
    glance(mod)$r.squared
  })) %>%
  dplyr::select(term, tidy, r2) %>%
  unnest()

# Temperature range and breeding range size r ~ 0.5

# write.csv(spp_traits, "model/spp_trait_model_output.csv", row.names = F)
spp_traits <- read.csv("model/spp_trait_model_output.csv", stringsAsFactors = F)

spp_traits_pred <- spp_traits %>%
  filter(p.value < 0.05)

### species traits models effect plots

trait_effects <- spp_table_traits %>%
  left_join(model_fits)

# propFor v deltaED
ggplot(filter(trait_effects, term == "deltaED"), aes(x = propFor, y = Estimate)) + 
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(y = "Effect of change in forest edge density", x = "Mean proportion forest cover")
ggsave("figures/main_analysis_figs/trait_model_deltaED_propFor.pdf")

# volume v tmin

ggplot(filter(trait_effects, term == "tmin"), aes(x = temp_range, y = Estimate)) + 
  geom_text(aes(label = SPEC)) + geom_smooth(method = "lm", se = F) +
  labs(y = "Effect of trend in Tmin", x = "Temperature range")
ggsave("figures/main_analysis_figs/trait_model_tmin_volume.pdf")

#### Range position models ####

## Mean breeding season temperature on BBS routes

route_climate <- read.csv("climate/bbs_routes_breeding_season_climate.csv", stringsAsFactors = F) %>%
  group_by(stateroute) %>%
  summarize(mean_tmax = mean(mean_tmax),
            mean_tmin = mean(mean_tmin))

## Join route position with clim_hab_poptrend table

clim_hab_poptrend_means <- clim_hab_poptrend_z %>%
  left_join(route_climate)

## Individual models of Tmin and Tmax with interactions with breeding season mean Tmin and Tmax over study period

 model_fits_position <- clim_hab_poptrend_means %>%
   group_by(aou, SPEC) %>%
   nest() %>%
   mutate(lmFit = map(data, ~{
     df <- .

     # get species weights matrix

     # Coordinates of BBS routes

     coords_bbs <- df %>%
       ungroup() %>%
       dplyr::select(stateroute) %>%
       distinct() %>%
       left_join(routes) %>%
       dplyr::select(stateroute, longitude, latitude)

     coords_mat <- as.matrix(coords_bbs[, -1])

     # Calculate nearest neighbors and make spatial neighborhood
     k0 <- knearneigh(coords_mat, longlat=T, k=1)
     k1 <- knn2nb(k0)

     # Find maximum neighbor distance and use this distance to define neighborhoods
     max.d <- max(unlist(nbdists(k1, coords_mat, longlat=T)))
     nb0.birds <- dnearneigh(coords_mat, 0, max.d, longlat=T)
     plot(nb0.birds, coords_mat)

     # Using a distance threshold of 100km
     nb1.birds <- dnearneigh(coords_mat,1,100,longlat=T)
     plot(nb1.birds, coords_mat)

     # Create spatial weights based on linear distance decay
     glist.birds <- nbdists(nb0.birds, coords_mat, longlat=T)
     glist.birds <- lapply(glist.birds, function(x) 1-(x/ceiling(max.d))) # Round up max.d so that all point get weights
     wt.birds <- nb2listw(nb0.birds, style='B', glist=glist.birds)

     spautolm(abundTrend ~ tmax*mean_tmax + tmin*mean_tmin, data = df, wt.birds, na.action = na.fail, family = "CAR")
   })) %>%
   mutate(nObs = map_dbl(data, ~{
     df <- .
     nrow(df)
   }))  %>%
   mutate(pred = purrr::map(lmFit, ~{
     mod <- .
     data.frame(pred.vals = mod$fit$fitted.values)
   })) %>%
   mutate(comb = map2(data, pred, ~bind_cols(.x, .y))) %>%
   mutate(r2 = map_dbl(comb, ~{
     df <- .
     summary(lm(abundTrend ~ pred.vals, df))$r.squared
   })) %>%
   mutate(lm_broom = map(lmFit,  ~{
     mod <- .
     sum <- summary(mod)
     df <- as.data.frame(sum$Coef)
     df$term = row.names(df)
     df
   })) %>%
   dplyr::select(aou, SPEC, nObs, lm_broom, r2) %>%
   unnest() %>%
   rename(std.error = `Std. Error`, z.value = `z value`, p.value = `Pr(>|z|)`)

# write.csv(model_fits_position, "model/range_position_model_tables.csv", row.names = F)
model_fits_position <- read.csv("model/range_position_model_tables.csv", stringsAsFactors = F)

model_fits_position_fig <- model_fits_position %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(p.value < 0.05 & p.value >= 0.01 ~ "0.01 < p < 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "0.001 < p < 0.01",
                         p.value < 0.001 ~ "p < 0.001",
                         TRUE ~ "p > 0.05")) %>%
  mutate(sig = fct_relevel(sig, c("p < 0.001", "0.001 < p < 0.01", "0.01 < p < 0.05", "p > 0.05"))) %>%
  mutate(dir = ifelse(Estimate > 0, "Positive effect", "Negative effect")) %>%
  ungroup() %>%
  filter(SPEC != "CERW")

# Figure

tmax_int_plot <- density_plot(model_fits_position_fig, "tmax:mean_tmax", "Interaction: trend in Tmax & mean Tmax")
tmin_int_plot <- density_plot(model_fits_position_fig, "tmin:mean_tmin", "Interaction: trend in Tmin & mean Tmin")
  
legend <- get_legend(tmax_int_plot)

temp_effects <- plot_grid(tmax_int_plot + theme(legend.position = "none"), 
                          tmin_int_plot + theme(legend.position = "none") + ylab(" "),
                          nrow = 1,
                          labels = c("A", "B"))
plot_grid(temp_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/range_position_temp_responses.pdf", units = "in", height = 3, width = 10)

## Supplemental figure: p-values

tmax_int_plot <- pval_plot(model_fits_position_fig, "tmax:mean_tmax", "Interaction: trend in Tmax & mean Tmax")
tmin_int_plot <- pval_plot(model_fits_position_fig, "tmin:mean_tmin", "Interaction: trend in Tmin & mean Tmin")

legend <- get_legend(tmax_int_plot)

temp_effects <- plot_grid(tmax_int_plot + theme(legend.position = "none"), 
                          tmin_int_plot + theme(legend.position = "none") + ylab(" "),
                          nrow = 1,
                          labels = c("A", "B"))
plot_grid(temp_effects, legend, rel_widths = c(2, 0.4))
ggsave("figures/main_analysis_figs/range_position_temp_pvals.pdf", units = "in", height = 3, width = 10)

#### Strongest species response for each predictor ####

# NOFL - deltaED

flicker <-  clim_hab_poptrend_z %>%
  filter(Common_name == "Northern flicker")

flicker_plot <- ggplot(flicker, aes(x = deltaED, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Change in edge density", y = "Abundance trend", title = "Northern flicker")

# INBU - deltaProp

bunting <- clim_hab_poptrend_z %>%
  filter(Common_name == "Indigo bunting") 

bunting_plot <- ggplot(bunting, aes(x = deltaProp, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Change in forest cover", y = "", title = "Indigo bunting")
 
# BAWW - tmax

warbler <- clim_hab_poptrend_z %>%
  filter(Common_name == "Black-and-white warbler")

warbler_plot <- ggplot(warbler, aes(x = tmax, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Trend in Tmax", y = "", title = "Black-and-white warbler")

# MAWA - tmin

magnolia <- clim_hab_poptrend_z %>%
  filter(Common_name == "Magnolia warbler")

magnolia_plot <- ggplot(magnolia, aes(x = tmin, y = abundTrend)) + geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Trend in Tmin", y = "Abundance trend", title = "Magnolia warbler")

# SUTA - deltaED:tmax

tanager <- clim_hab_poptrend_z %>%
  filter(Common_name == "Summer tanager") 

tanager$tmax_sign <- ifelse(tanager$tmax < 0, "Weakest warming trends", "Strongest warming trends")
tanager_plot <- ggplot(tanager, aes(x = deltaED, y = abundTrend)) + geom_point() + 
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~tmax_sign) +
  theme(panel.spacing = unit(4, "lines")) +
  labs(x = "Change in edge density", y = "Abundance trend", title = "Summer tanager")

## Figure for MS

indiv_spp_add <- plot_grid(flicker_plot, bunting_plot, magnolia_plot, warbler_plot, nrow = 2, labels = c("A", "B", "C", "D"))
indiv_spp_multi <- plot_grid(indiv_spp_add, tanager_plot, nrow = 2, rel_heights = c(0.66, 0.33),
                             labels = c(" ", "E", "F"))
ggsave("figures/main_analysis_figs/indiv_spp_multipanel.pdf", indiv_spp_multi, units = "in", width = 9, height = 10)

