

# list packages that will be needed
pacakges <- c("sp","raster","rgeos","maptools","rworldmap", "rgdal","maps", 
              "factoextra","RStoolbox", "corrplot", "fpc","RColorBrewer", "Rfast",
              "rpart", "rpart.plot" ,"NbClust",  "psych"  , "reshape2", "tidyverse", "magrittr")

# install those packages that you dont have yet
install.packages( 
  setdiff(pacakges,installed.packages()[,1] )
)




# clear environment
rm(list = ls())

gc()



#Load packages:

library(raster)
library(rgeos)
library(maptools)
library(factoextra)
library(RStoolbox)
library(reshape2)
library(corrplot)
library(fpc)
library(RColorBrewer)
library(rpart)
library(rpart.plot)
library(NbClust)
library(sp)
library(psych)
library(rworldmap)
library(Rfast)
library(rgdal)
library(tidyverse)
library(magrittr)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                                 __Load data__                                            ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 




wd <-"C:/Users/Nils/Documents/R/MBERC - Practicals and projects/Kirsty"
setwd(wd)



paste0(wd,"/inputs") -> rasterDir
# Add your own rasters in there

paste0(wd,"/outputs")  -> outputsDir
# in case you dont have it yet, create the directory
outputsDir %>%  dir.create()

paste0(outputsDir,"/NA") -> NADir
# in case you dont have it yet, create the directory
NADir %>%  dir.create()

paste0(outputsDir,"/HA") -> HADir
HADir %>%  dir.create()


#Load environmental variables: # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# make a raster names list
rasters <- list.files(path =  rasterDir,
                      pattern = ".tif$",
                      full.names = TRUE) 
#make a vector of rasters names
rasters_names <- list.files(path =  rasterDir,
                            pattern = ".tif$",
                            full.names = FALSE)  

# use one of the larers extent to set the stack's extent
ext <- extent(raster(str_subset(rasters,pattern = "2006-2015_fluxXYZ112deg.tif"))) 


# make a raster stack
var <- stack()

for (i in seq(rasters)) {
  # load rasters and crop to same extent
  addLayer(var, raster(rasters[i]) %>%  crop(ext))  -> var
}# next raster

# give shorter names to individual layers
namestable <- tibble(
  rasters_files = rasters_names,
  names = c(
    "poc",
    "salinity",
    "temperature",
    "bbpi",
    "depth",
    "fbpi",
    "slope",
    "currvel",
    "dissox",
    "nit",
    "phos",
    "sil"
  ), 
  chicNames = c(
    'poc' ,
    'Salinity' ,
    'Temperature' ,
    'BBPI' ,
    'Depth',
    'FBPI',
    'Slope',
    'Current velocity',
    'Dissolved oxygen',
    'Nitrate',
    'Phosphate',
    'Silicate' 
  )
)

names(var) <- namestable$names


# make a mask of oceans (elevation below sea-level )
var$depth -> depth # use the depth layer
var %<>% dropLayer("depth") # remove depth from stack

mask_ <- depth <= 0
mask <- crop(mask_, var)
plot(mask)



# load world continents for nicer maps
rworldmap::getMap(resolution = "coarse") -> continents

# load the rasters of the clustering
paste0(outputsDir,"/","Level1.tif") %>% raster() -> clara.rast
names(clara.rast) <- 'zone'


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                                 __Confusion Index__                                      ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 


#calculating confusion index for Level 1 clustering of the hierarchical approach. 
#Confusion index values range from 0 to 1, with values approaching 1 indicating greater 
#uncertainty in the clara clustering solution between two or more clusters.

#Convert rasters to vectors: - - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - - 


data <- tibble("ID" = 1:ncell(var)) 

vectors <- list()
for (v in var@layers) {
  v %>% 
    as.matrix() %>%
    as.vector() %>% 
    bind_cols(data, .) -> data
}




# turn the 0s into NAs
data[data == 0] <- NA
data_na <- na.omit(data)
# remove the ID column
data_norm <-
  (scale(select(data_na,-ID), scale = TRUE)) 


#Run Principal Component Analysis and clustering:

# PCA:
data_PC <- prcomp(data_norm, scale = FALSE)
eigen <- data_PC$sdev ^ 2
eigen # Examine eigenvalues to identify number of PCs with eigenvalues > 1
dataPC_rotated <-
  psych::principal(data_norm,
                   rotate = "varimax",
                   nfactors = 4,
                   scores = TRUE) # Change nfactors to desired number of PCs
data_PCs <- dataPC_rotated$scores

# Clustering:
fit <- kmeans(data_PCs, 4, iter.max = 40)

# Create matrix of input data & clustering output:

m <- matrix(NA, nrow(data_PCs), ncol = nrow(fit$centers))

for (i in 1:nrow(data_PCs)) {
  m[i, ] <- as.matrix(Dist(rbind(data_PCs[i, ],  fit$centers)))[-1, 1]
}

mu <- matrix(NA, nrow(data_PCs), ncol = nrow(fit$centers))

# Calculate confusion index:
for (i in 1:nrow(data_PCs)) {
  mu[i, ] <- 1 / m[i, ] ^ 2 * 1 / sum(1 / m[i, ] ^ 2)
}

CI <- numeric(nrow(data_PCs))

for (i in 1:nrow(data_PCs)) {
  CI[i] <-
    mu[i, order(mu[i, ], decreasing = TRUE)[2]] / mu[i, order(mu[i, ], decreasing =
                                                                TRUE)[1]]
}

data_PCs2 <- data.frame(data_PCs, CI)
temp1 <- cbind(data_na, data_PCs2$CI)



# a column for confidence indices
temp1 %>%  as.data.frame() %>%  as_tibble() %>%
  select(ID, `data_PCs2$CI`) -> d

# remake the raster and include the NA cells that had been removed
matrix1 <- matrix(nrow = nrow(data), ncol = 2)
# join the CI values with the cells that have one
matrix1 %>% as_tibble() %>%
  mutate(ID = 1:nrow(.)) %>%
  left_join(d, by = "ID") %>%
  rename(CI = `data_PCs2$CI`) -> dd

# make raster at the dimension of the input raster
matrix3 <- matrix(dd$CI, nrow = nrow(var), ncol = ncol(var))

#Plot and save output as raster:
Level1_CI <- raster(
  matrix3,
  xmn = -180,
  xmx = 180,
  ymn = -90,
  ymx = 90
)
plot(Level1_CI)
crs(Level1_CI) <- "+init=epsg:4326"
writeRaster(Level1_CI, filename = paste0(HADir, "/" , "Level1_CI.tif"))

System.time()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                         #__Crop outputs to South Atlantic study area:__                                  ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



#Load mask of South Atlantic:
SAmask <-
  readOGR(paste0(rasterDir , "/", "Mask_SouthAtlantic.shp"))
#Crop classification Level 1: - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -

level1 <- clara.rast
level1.crop <- crop(level1, SAmask)
level1 <- mask(level1.crop, SAmask)
# plot map of the results results
plot(level1)
plot(continents, add = T, col = "grey")
# Export
writeRaster(level1, paste0(HADir, "/", "Level1_SAtlant.tif"), overwrite =
              TRUE)

#Crop confidence map Level 1: - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -
Level1_CI.crop <- crop(Level1_CI, SAmask)
Level1_CI <- mask(Level1_CI.crop, SAmask)

# plot map of the results results
plot(Level1_CI)
plot(continents, add = T, col = "grey")
# Export
writeRaster(Level1_CI,
            paste0(HADir, "/", "Level1_CI_SAtlant.tif"),
            overwrite = TRUE)



