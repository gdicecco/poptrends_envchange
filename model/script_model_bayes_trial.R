### Hierarchical Bayes approach to modeling abundance trends in response to climate and land use change

library(tidyverse)
library(rstan)
library(brms)

options(mc.cores = 4)

## Data to use
bbs_env <- read_csv("model/bbs_route_env_change.csv")

spp_table_traits <- read.csv("traits/forest_spp_traits_MS.csv", stringsAsFactors = F)

bbs_subset <- read_csv("model/BBS_counts_model.csv") %>%
  left_join(dplyr::select(bbs_env, -year), by = c("stateroute")) %>%
  na.omit() %>%
  filter(aou %in% spp_table_traits$AOU) %>%
  arrange(aou, stateroute, year) %>%
  group_by(aou, stateroute) %>%
  mutate(n_year = n_distinct(year)) %>%
  filter(n_year >= 10) %>%
  mutate(st_aou = paste0(aou, stateroute),
    aou_index = as.numeric(as.factor(aou)),
    index = as.numeric(as.factor(st_aou))) ## I()


## BRMS test level 1 abundance trends
# fit <- brm(speciestotal ~ year + first_yr + (year + first_yr|st_aou),
#            data = bbs_subset,
#            family = poisson,
#            inits = "0") 
# 
# make_stancode(speciestotal ~ year + first_yr + (year + first_yr|st_aou),
#               data = bbs_subset,
#               family = poisson)

## Setup stan data
Y_vec <- bbs_subset$speciestotal ## Counts vec to predict
n_spp <- length(unique(bbs_subset$aou)) ## Number of spp
n_obs <- length(Y_vec) ## Total number of points
n_spp_sites <- length(unique(bbs_subset$st_aou))

spp_vec <- bbs_subset$aou_index
year_vec <- bbs_subset$year ## Years

site_data <- bbs_subset %>%
  ungroup() %>%
  select(st_aou, tmax, tmin, deltaED) %>%
  distinct()

edge_vec <- site_data$deltaED
tmax_vec <- site_data$tmax
tmin_vec <- site_data$tmin

nst_vec <- bbs_subset$index ## Species x site indexing

st_aous <- unique(bbs_subset$st_aou)
nsp_vec <- as.numeric(as.factor(substr(st_aous, 1, 4))) ## species IDs for species/site index

## trait data
trait_data <- spp_table_traits %>%
  arrange(AOU)

for_vec <- trait_data$propFor
temp_vec <- trait_data$temp_range
area_vec <- trait_data$Brange_Area_km2

## List all data
abund_data <- list(y = Y_vec, 
                   Nsp = n_spp, 
                   Nst = n_spp_sites,
                   N = n_obs,
                   sp = spp_vec,
                   cn_id = nst_vec,
                   cn_sp = nsp_vec,
                   year = year_vec,
                   edge = edge_vec,
                   tmax = tmax_vec,
                   tmin = tmin_vec,
                   propFor = for_vec,
                   temp_range = temp_vec,
                   Brange_Area_km2 = area_vec)

n_sam <- 10000
n_warmup <- round(n_sam/5,0)
n_chain <- 4

abund_fit <- stan(file='model/bbs_abund_mod_env.stan',
                  data = abund_data,
                  iter =(n_sam+n_warmup), 
                  warmup=n_warmup,
                  cores = 4,
                chains = n_chain,
              control = list(adapt_delta = 0.96, max_treedepth = 12, stepsize = 0.004))

saveRDS(abund_fit, "model/abund_bayes_fit.rds")

 