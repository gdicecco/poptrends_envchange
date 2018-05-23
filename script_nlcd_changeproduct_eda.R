library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("\\\\BioArk\\hurlbertlab\\GIS\\LandCoverData\\")

setwd("\\\\BioArk\\hurlbertlab\\GIS\\LandCoverData\\nlcd_landcover_change\\")
codes <- read.csv("nlcd_1992_to_2001_landcover_change\\anderson_land_cover_codes.csv", stringsAsFactors = F) # NLCD land cover class codes

## Combine NLCD regions into one for 1992 - 2001
setwd("nlcd_1992_to_2001_landcover_change")
files <- list.files()
area.files <- files[str_detect(files, "area")]
dir <- getwd()

# Get area changeproduct raster file from NLCD directory
get.file <- function(x) {
  files2 <- list.files(paste0(dir, "/", area.files[x], ""))
  file.path <- files2[str_detect(files2, "area")]
  return(list(folder = area.files[x], file.name = file.path))
}

# Function to merge two raster files with different origins and extents, same resolution and crs
merge.areas <- function(x, y) {
  extent1 <- extent(x)
  extent2 <- extent(y)
  z <- raster(xmn = min(extent1[1], extent2[1]), xmx = max(extent1[2], extent2[2]), # make blank raster file with extent of the two files added
              ymn = min(extent1[3], extent2[3]), ymx = max(extent1[4], extent2[4]), 
              resolution = c(900, 900), crs = crs(x))
  merge1 <- mosaic(z, x, fun = max, tolerance = max(abs(origin(z) - origin(x)))) # function: where rasters overlap keep larger value
  merge2 <- mosaic(merge1, y, fun = max, tolerance = max(abs(origin(merge1) - origin(y)))) # tolerance allows for different origins 
  return(merge2) # returns the two rasters merged into one
}

# Loop to merge changeproduct areas

file1 <- get.file(1)
data <- raster(paste0(dir, "/", file1$folder, "/", file1$file.name, sep = ""))
data.km <- aggregate(data, fact = 30, fun = modal) #modal - mode of the raster values in the cell

file2 <- get.file(2)
data2 <- raster(paste0(dir, "/", file2$folder, "/", file2$file.name, sep = ""))
data2.km <- aggregate(data2, fact = 30, fun = modal) #modal - mode of the raster values in the cell

region <- merge.areas(data.km, data2.km)
crs <- crs(region)

for(i in 3:length(area.files)) { 
  if(i == 8) { # area 8 - Michigan - has different projection
    file3 <- get.file(i)
    data3 <- raster(paste0(dir, "/", file3$folder, "/", file3$file.name, sep = ""))
    data3.km <- aggregate(data3, fact = 30, fun = modal) #modal - mode of the raster values in the cell
    data3.proj <- projectRaster(data3.km, crs = crs)
    region <- merge.areas(region, data3.proj)
  } else {
  file3 <- get.file(i)
  data3 <- raster(paste0(dir, "/", file3$folder, "/", file3$file.name, sep = ""))
  data3.km <- aggregate(data3, fact = 30, fun = modal) #modal - mode of the raster values in the cell
  region <- merge.areas(region, data3.km)
  }
}
plot(region)

## Read in combined 1992-2001 US raster
region <- raster("1992-2001_changeproduct_US.grd")

## Stack with BCRs
# read in BCR shape file, convert to raster

bcrshp <- readOGR("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\bcr_terrestrial_shape\\BCR_Terrestrial_master.shp") #BCRs
bcr.naomit <- bcrshp[!is.na(bcrshp@data$REGION), ] # remove NAs
usa <- bcr.naomit[bcr.naomit@data$REGION == "USA", ] # USA only
contig.us <- usa[!(usa@data$PROVINCE_S == "ALASKA" | usa@data$PROVINCE_S == "HAWAIIAN ISLANDS"), ] # continental US only
bcrs.us <- contig.us@data$BCR

# need to re-project contig.us to match crs of region
us.proj <- sp::spTransform(contig.us, crs(region))

blank <- raster(ext = extent(region), crs = crs(region), resolution = res(region)) # empty raster
raster.us <- rasterize(us.proj, blank, field = "BCR", fun = "first") # raster of BCRs

# stack BCRs with NLCD data
bcr.nlcd <- stack(raster.us, region)
names(bcr.nlcd) <- c("BCRs", "NLCD")
plot(bcr.nlcd)

# identify codes that are categories of interest 
# (natural --> urban (forest - 42, grassland - 52), natural --> agriculture (forest - 46, grassland - 56))
# compare across BCR regions to find areas with high amounts of these at coarse scale 

# Raster to data frame
stack.df <- as.data.frame(rasterToPoints(bcr.nlcd))

summary <- stack.df %>%
  filter(NLCD == 42 | NLCD == 52 | NLCD == 46 | NLCD == 56) %>%
  group_by(BCRs, NLCD) %>%
  summarize(total = n())

theme_set(theme_bw())
code.names <- data.frame(nlcd = c(42, 52, 46, 56),
                         change = c("Forest-Urban", "Grassland-Urban","Forest-Agriculture", "Grassland-Agriculture"))
summary.plot <- left_join(summary, code.names, by = c("NLCD" = "nlcd"))
plot1992 <- ggplot(summary.plot, aes(x = BCRs, y = total, fill = BCRs)) + geom_col() + facet_wrap(~change) + ylab("No. 900 m x 900 m cells") + ggtitle("1992-2001")
# 1992-2001 land cover changes of interest per BCR

# 2001-2006
setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\GIS\\LandCoverData\\nlcd_landcover_change\\")
files <- list.files()
nlcd.files <- files[str_detect(files, "2006")]
dir <- getwd()

get.file.img <- function(x) {
  files2 <- list.files(paste0(dir, "/", nlcd.files[x], ""))
  file.path <- files2[str_detect(files2, "img")]
  return(list(folder = nlcd.files[x], file.name = file.path))
}

file.2001 <- get.file.img(1)
nlcd2001 <- raster(paste0(dir, "/", file.2001$folder, "/", file.2001$file.name, sep = ""))

nlcd2001.df <- as.data.frame(levels(nlcd2001))
toAnthro <- as.data.frame(levels(nlcd2001)) %>%
  filter(grepl("Developed|Pasture|Crops", X2006.Class)) %>%
  filter(grepl("Forest|Grassland|Shrub", X2001.Class))
# Filtered change pixels to pixels that went from natural areas (forest, grassland, shrubland) to used by humans (urban, agriculture)

nlcd2001.km <- aggregate(nlcd2001, fact = 30, fun = modal)

us.proj <- sp::spTransform(contig.us, crs(nlcd2001.km))

blank <- raster(ext = extent(nlcd2001.km), crs = crs(nlcd2001.km), resolution = res(nlcd2001.km)) # empty raster
raster.us <- rasterize(us.proj, blank, field = "BCR", fun = "first") # raster of BCRs

# stack BCRs with NLCD data
bcr.nlcd.2001 <- stack(raster.us, nlcd2001.km)
names(bcr.nlcd.2001) <- c("BCRs", "NLCD")
plot(bcr.nlcd.2001)

# Raster to data frame
stack.df <- as.data.frame(rasterToPoints(bcr.nlcd.2001))

# group pixels four categories: forest -> ag, forest -> dev, grass/shrub -> ag, grass/shrub -> dev
code.names <- data.frame(change = c(13, 23, 14, 24),
                         class = c("Forest-Urban", "Grassland-Urban","Forest-Agriculture", "Grassland-Agriculture"))

summary01 <- stack.df %>%
  filter(NLCD %in% toAnthro$ID) %>%
  left_join(toAnthro, by = c("NLCD" = "ID")) %>%
  mutate(from = ifelse(grepl("Forest", X2001.Class), 1, 2),
         to = ifelse(grepl("Developed", X2006.Class), 3, 4),
         change = from*10+to) %>% 
  group_by(BCRs, change) %>%
  summarize(total = n()) %>%
  left_join(code.names, by = "change")

theme_set(theme_bw())
plot2001 <- ggplot(summary01, aes(x = BCRs, y = total, fill = BCRs)) + geom_col() + facet_wrap(~class) + ylab("No. 900 m x 900 m cells") + ggtitle("2001-2006")

# 2006-2011
file.2006 <- get.file.img(2)
nlcd2006 <- raster(paste0(dir, "/", file.2006$folder, "/", file.2006$file.name, sep = ""))

nlcd2006.df <- as.data.frame(levels(nlcd2006))
toAnthro <- as.data.frame(levels(nlcd2006)) %>%
  filter(grepl("Developed|Pasture|Crops", X2011.Class)) %>%
  filter(grepl("Forest|Grassland|Shrub", X2006.Class))
# Filtered change pixels to pixels that went from natural areas (forest, grassland, shrubland) to used by humans (urban, agriculture)

nlcd2006.km <- aggregate(nlcd2006, fact = 30, fun = modal)

bcr.nlcd.2006 <- stack(raster.us, nlcd2006.km)
names(bcr.nlcd.2006) <- c("BCRs", "NLCD")
plot(bcr.nlcd.2006)

# Raster to data frame
stack.2006 <- as.data.frame(rasterToPoints(bcr.nlcd.2006))

summary2006 <- stack.2006 %>%
  filter(NLCD %in% toAnthro$ID) %>%
  left_join(toAnthro, by = c("NLCD" = "ID")) %>%
  mutate(from = ifelse(grepl("Forest", X2006.Class), 1, 2),
         to = ifelse(grepl("Developed", X2011.Class), 3, 4),
         change = from*10+to) %>% 
  group_by(BCRs, change) %>%
  summarize(total = n()) %>%
  left_join(code.names, by = "change")

theme_set(theme_bw())
plot2006 <- ggplot(summary2006, aes(x = BCRs, y = total, fill = BCRs)) + geom_col() + facet_wrap(~class) + ylab("No. 900 m x 900 m cells") + ggtitle("2006-2011")

# Land cover change over three time windows for each BCR

library(cowplot)
multiplot <- plot_grid(plot1992, plot2001, plot2006, nrow = 2)

# How much conversion of natural to human use land in each strata for three time windows
library(tidyr)
colnames(summary.plot)[c(2,4)] <- c("change", "class")
totalchange <- rbind(data.frame(year = 1992, summary.plot),
                     data.frame(year = 2001, summary01),
                     data.frame(year = 2006, summary2006)) %>%
  group_by(year, BCRs) %>%
  summarize(total.fragment = sum(total)) %>%
  group_by(year) %>%
  nest()

# Subset BCR shapefile (us.proj) to just BCR column, use amount of conversion to anthropogenic land use to color strata

us.bcrs <- us.proj[, -c(2:8)]
bcrs <- unique(us.bcrs@data$BCR)
bcrs.df <- data.frame(BCRs = bcrs, total.fragment = 0)

library(purrr)
bcrs.frag <- map(totalchange$data, ~{
  right_join(., bcrs.df)
})

bcrs.frag <- map(totalchange$data, ~{
  left_join(us.bcrs@data, ., by = c("BCR" = "BCRs"))
})

us.bcrs@data$Frag92 <- bcrs.frag[[1]]$total.fragment
us.bcrs@data$Frag01 <- bcrs.frag[[2]]$total.fragment
us.bcrs@data$Frag06 <- bcrs.frag[[3]]$total.fragment

us.bcrs@data$id <- rownames(us.bcrs@data)
bcr.points <- fortify(us.bcrs, region="id")
bcr.df <- plyr::join(bcr.points, us.bcrs@data, by="id")

# Plots with just land cover change 

map92 <- ggplot(bcr.df) + aes(long, lat, group = group, fill = Frag92) + geom_polygon() + geom_path(color = "black") + coord_equal() +
  scale_fill_continuous(low = "white", high = "red", limits = c(1, 3000)) + ggtitle("1992-2001") + theme(axis.line = element_blank(),
                                                                                    axis.title = element_blank(),
                                                                                    axis.text = element_blank(),
                                                                                    axis.ticks = element_blank()) + labs(fill = "No. pixels")

map01 <- ggplot(bcr.df) + aes(long, lat, group = group, fill = Frag01) + geom_polygon() + geom_path(color = "black") + coord_equal() +
  scale_fill_continuous(low = "white", high = "red", limits = c(1, 3000)) + ggtitle("2001-2006") + theme(axis.line = element_blank(),
                                                                                    axis.title = element_blank(),
                                                                                    axis.text = element_blank(),
                                                                                    axis.ticks = element_blank()) + labs(fill = "No. pixels")

map06 <- ggplot(bcr.df) + aes(long, lat, group = group, fill = Frag06) + geom_polygon() + geom_path(color = "black") + coord_equal() +
  scale_fill_continuous(low = "white", high = "red", limits = c(1, 3000)) + ggtitle("2006-2011") + theme(axis.line = element_blank(),
                                                                                    axis.title = element_blank(),
                                                                                    axis.text = element_blank(),
                                                                                    axis.ticks = element_blank()) + labs(fill = "No. pixels")
plot_grid(map92, map01, map06, nrow = 2)

# Knit together BBS route paths and overlay onto data about land cover transitions

setwd("\\\\Bioark.bio.unc.edu/hurlbertlab/Databases/BBS/GPS_stoplocations/")

us_routes <- readOGR("bbsrte_2012_alb/bbsrte_2012_alb.shp")

# subset routes that are between 38000 and 42000 m, remove Alaska (rteno between 3000 and 4000)
# reproject to crs of us.bcrs
# use fortify to convert to df

us_routes_short <- us_routes[us_routes@data$rte_length < 42000 & us_routes@data$rte_length > 38000, ]
us_subs <- us_routes_short[!(us_routes_short@data$rteno < 4000 & us_routes_short@data > 3000), ]

# crs.us <- CRS(us.bcrs)
us_subs <- spTransform(us_subs, crs.us)

us_subs@data$id <- rownames(us_subs@data)
us_routes.df <- fortify(us_subs, region = "id")

# Plots with routes 

routes92 <- ggplot() + geom_polygon(bcr.df, aes(long, lat, group = group, fill = Frag92)) + coord_equal() +
  geom_path(bcr.df, aes(long, lat, group = group), color = "gray19") +
  geom_path(us_routes.df, aes(long, lat, group = group), color = "black") +
  scale_fill_continuous(low = "white", high = "red", limits = c(1, 3000)) + ggtitle("1992-2001") + theme(axis.line = element_blank(),
                                                                                                         axis.title = element_blank(),
                                                                                                         axis.text = element_blank(),
                                                                                                         axis.ticks = element_blank()) + labs(fill = "No. pixels")

routes01 <- ggplot() + geom_polygon(bcr.df, aes(long, lat, group = group, fill = Frag01)) + coord_equal() +
  geom_path(bcr.df, aes(long, lat, group = group), color = "gray19") +
  geom_path(us_routes.df, aes(long, lat, group = group), color = "black") +
  scale_fill_continuous(low = "white", high = "red", limits = c(1, 3000)) + ggtitle("2001-2006") + theme(axis.line = element_blank(),
                                                                                                         axis.title = element_blank(),
                                                                                                         axis.text = element_blank(),
                                                                                                         axis.ticks = element_blank()) + labs(fill = "No. pixels")

routes06 <- ggplot() + geom_polygon(bcr.df, aes(long, lat, group = group, fill = Frag06)) + coord_equal() +
  geom_path(bcr.df, aes(long, lat, group = group), color = "gray19") +
  geom_path(us_routes.df, aes(long, lat, group = group), color = "black") +
  scale_fill_continuous(low = "white", high = "red", limits = c(1, 3000)) + ggtitle("2006-2011") + theme(axis.line = element_blank(),
                                                                                                         axis.title = element_blank(),
                                                                                                         axis.text = element_blank(),
                                                                                                         axis.ticks = element_blank()) + labs(fill = "No. pixels")



