## NLCD Changeproduct explore plots
library(tidyverse)

setwd("C:/Users/gdicecco/Desktop/data/nlcd_fragmentation_longleaf_out/")
files <- list.files()

# Create master data frame of all individual fragmenation indices csvs
nlcd <- data.frame(file = NA, year = NA, bcr = NA)
for(i in 1:length(files)) {
  if(i == 1) {
    name <- files[1]
    data <- read.csv(name, stringsAsFactors = F)
    data$class <- round(data$class)
    data.short <- filter(data, class > 0)
    nlcd[i:nrow(data.short), "file"] <- name
    nlcd[i:nrow(data.short), "year"] <- as.numeric(str_match_all(name, "[0-9]+")[[1]][3,1])
    nlcd[i:nrow(data.short), "bcr"] <- as.numeric(str_match_all(name, "[0-9]+")[[1]][4,1])
    nlcd <- bind_cols(nlcd, data.short)
  } else {
    name <- files[i]
    data <- read.csv(name, stringsAsFactors = F)
    data$class <- round(data$class)
    data.short <- filter(data, class > 0)
    data.short$file <- name
    data.short$year <- as.numeric(str_match_all(name, "[0-9]+")[[1]][3,1])
    data.short$bcr <- as.numeric(str_match_all(name, "[0-9]+")[[1]][4,1])
    nlcd <- bind_rows(nlcd, data.short)
  }
}
setwd("//BioArk/HurlbertLab/DiCecco/data/")
write.csv(nlcd, "fragmentation_indices_nlcd.csv", row.names = F)
