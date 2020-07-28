# Create master dataset from route level fragmentation output

library(tidyverse)

setwd("C:/Users/gdicecco/Desktop/data/canada_frag_output/")
files <- list.files()

# Create master data frame of all individual fragmenation indices csvs
nlcd <- data.frame(file = NA, year = NA, stateroute = NA)
for(i in 1:length(files)) {
  if(i == 1) {
    name <- files[i]
    data <- read.csv(name, stringsAsFactors = F)
    data$class <- round(data$class)
    data.short <- filter(data, class > 0)
    nlcd[i:nrow(data.short), "file"] <- name
    nlcd[i:nrow(data.short), "year"] <- as.numeric(str_match_all(name, "[0-9]+")[[1]][3,1])
    nlcd[i:nrow(data.short), "stateroute"] <- as.numeric(str_match_all(name, "[0-9]+")[[1]][6,1])
    nlcd <- bind_cols(nlcd, data.short)
  } else {
    name <- files[i]
    data <- read.csv(name, stringsAsFactors = F)
    data$class <- round(data$class)
    if(unique(data$class) != 0){
    data.short <- filter(data, class > 0)
    data.short$file <- name
    data.short$year <- as.numeric(str_match_all(name, "[0-9]+")[[1]][3,1])
    data.short$stateroute <- as.numeric(str_match_all(name, "[0-9]+")[[1]][6,1])
    nlcd <- bind_rows(nlcd, data.short)
    } else {print(name)}
    
  }
}

setwd("\\\\BioArk\\hurlbertlab\\DiCecco\\data\\")
write.csv(na.omit(nlcd), "fragmentation_indices_canada.csv", row.names = F)
