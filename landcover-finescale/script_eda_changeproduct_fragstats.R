## NLCD Changeproduct explore plots
library(tidyverse)

setwd("C:/Users/gdicecco/Desktop/csvs/")
files <- list.files()
files <- files[!grepl("routes_routes", files)]

# Create master data frame of all individual fragmenation indices csvs
nlcd <- data.frame(file = NA, year1 = NA, year2 = NA, bcr = NA, extent = NA)
for(i in 1:length(files)) {
  if(i == 1) {
  name <- files[i]
  data <- read.csv(name, stringsAsFactors = F)
  data$class <- round(data$class)
  data.short <- filter(data, class > 0)
  nlcd[i:nrow(data.short), "file"] <- name
  nlcd[i:nrow(data.short), "year1"] <- as.numeric(str_match_all(name, "[0-9]+")[[1]][3,1])
  nlcd[i:nrow(data.short), "year2"] <- as.numeric(str_match_all(name, "[0-9]+")[[1]][4,1])
  nlcd[i:nrow(data.short), "bcr"] <- as.numeric(str_match_all(name, "[0-9]+")[[1]][5,1])
  nlcd[i:nrow(data.short), "extent"] <- ifelse(grepl(name, "routes"), "routes", "bcr")
  nlcd <- bind_cols(nlcd, data.short)
  } else {
    name <- files[i]
    data <- read.csv(name, stringsAsFactors = F)
    data$class <- round(data$class)
    data.short <- filter(data, class > 0)
    data.short$file <- name
    data.short$year1 <- as.numeric(str_match_all(name, "[0-9]+")[[1]][3,1])
    data.short$year2 <- as.numeric(str_match_all(name, "[0-9]+")[[1]][4,1])
    data.short$bcr <- as.numeric(str_match_all(name, "[0-9]+")[[1]][5,1])
    data.short$extent <- ifelse(grepl(name, "routes"), "routes", "bcr")
    nlcd <- bind_rows(nlcd, data.short)
  }
}
setwd("//BioArk/HurlbertLab/DiCecco/data/")
write.csv(nlcd, "fragmentation_indices_nlcd_changeproducts.csv", row.names = F)

##### EDA of fragmentation indices ####
setwd("C:/Users/gdicecco/Desktop/git/NLCD_fragmentation/landcover-finescale/")

