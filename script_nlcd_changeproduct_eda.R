library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("C:/Users/gdicecco/Desktop/1992_2001_nlcd/")
codes <- read.csv("anderson_land_cover_codes.csv", stringsAsFactors = F) # NLCD land cover class codes

# Combine NLCD regions into one for 1992 - 2001
setwd("C:/Users/gdicecco/Desktop/1992_2001_nlcd/")
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

setwd("C:/Users/gdicecco/Desktop/1992_2001_nlcd/")
writeRaster(region, "1992-2001_changeproduct_US_noMI.grd", overwrite = T)

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
ggplot(summary.plot, aes(x = BCRs, y = total, fill = BCRs)) + geom_col() + facet_wrap(~change) + ylab("No. 900 m x 900 m cells")
# 1992-2001 land cover changes of interest per BCR

# Other time points (2001-2006, 2006-2011)
setwd("\\\\Bioark.bio.unc.edu\\hurlbertlab\\GIS\\LandCoverData\\nlcd_landcover_change\\")
files <- list.files()
nlcd.files <- files[str_detect(files, "2006")]
dir <- getwd()

get.file.img <- function(x) {
  files2 <- list.files(paste0(dir, "/", nlcd.files[x], ""))
  file.path <- files2[str_detect(files2, "img")]
  return(list(folder = area.files[x], file.name = file.path))
}

file.2001 <- get.file.img(1)
nlcd2001 <- raster(paste0(dir, "/", file.2001$folder, "/", file.2001$file.name, sep = ""))
nlcd2001.km <- aggregate(nlcd2001, fact = 30, fun = modal)
nlcd2001.df <- as.data.frame(rasterToPoints(nlcd2001.km)) # Need to find out what these indexes mean - 0 to 289

# Land cover change over three time windows for each BCR

# Eventually: knit together BBS route paths and overlay onto data about land cover transitions
