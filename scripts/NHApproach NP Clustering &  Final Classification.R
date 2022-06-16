title: "Broad-scale benthic habitat classification of the South Atlantic"
author: "Tiago Gandra, Kirsty McQuaid, Oli Hogg & Nils Piechaud"
date: "April 2022"


# =================================================================================================================
##                                  Non Hierarchichal clustering Approaching                                 ======
# =================================================================================================================

# list packages that will be needed
pacakges <- c("sp","raster","rgeos","maptools","rworldmap", "rgdal","maps", 
              "factoextra","RStoolbox", "corrplot", "fpc","RColorBrewer", "Rfast",
              "rpart", "rpart.plot" ,"NbClust",  "psych"  , "reshape2", "tidyverse", "magrittr")

# install those packages that you dont have yet
install.packages( 
  setdiff(pacakges,installed.packages()[,1] )
)

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

# enter the pathway
wd <-"C:/your/pathway/here"
# if you are using the Github repo in an Rstudio project, the WD is set for you
# wd <- getwd() # get the path to the the WD

setwd(wd)


paste0(wd,"/inputs") -> rasterDir
# Add your own rasters in there

paste0(wd,"/outputs")  -> outputsDir
# in case you dont have it yet, create the directory
outputsDir %>%  dir.create()

paste0(outputsDir,"/NA") -> NADir
# in case you dont have it yet, create the directory
NADir %>%  dir.create()

paste0(outputsDir,"/NHA") -> NHADir
# in case you dont have it yet, create the directory
NHADir %>%  dir.create()

#Load environmental variables: # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

# load world continents for nicer maps
rworldmap::getMap(resolution = "coarse") -> continents

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                                 pre processing                                          ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# subset the raster stack 
var %<>%  raster::subset(c("bbpi","fbpi", "slope", "temperature", "salinity", "poc")) 

# Rescale variables to have equal variance and a common scale of 0-1:

cor_Rasters <- stack()

for (i in names(var) %>%  seq()) {
  
  print(names(var)[i])
  
  r <- var[[i]]
  
  mean <- (cellStats(r, stat='mean'))
  stdev <- (cellStats(r, stat='sd'))
  norm <- ((r-mean)/(stdev))
  min <- (cellStats(norm, stat='min'))
  max <- (cellStats(norm, stat='max'))
  cor <- ((norm-min)/(max-min))
  
  cor_Rasters %<>% addLayer(cor)
}

# Extract to points, combine, remove NAs and georeferencing: = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  
# Topography:
bDF <- as_tibble(rasterToPoints(cor_Rasters[["bbpi"]]))
fDF <- as_tibble(rasterToPoints(cor_Rasters[["fbpi"]]))
slDF <- as_tibble(rasterToPoints(cor_Rasters[["slope"]]))

# merge all topography layers by their XY coordinates
topo_DF <- left_join(bDF, fDF, by = c("x", "y")) %>%
  left_join(slDF, by = c("x", "y"))

topo_XY <- topo_DF %>% na.omit() %>%  select(x, y)
topo_DF %<>%  na.omit(tibble(.$bbpi, .$fbpi, .$slope)) %>% select(-x,-y) # Add all variables to new data frame and select only columns with data in to remove georeferencing

# Water mass structure:
tDF <- as_tibble(rasterToPoints(cor_Rasters[["temperature"]]))
saDF <- as_tibble(rasterToPoints(cor_Rasters[["salinity"]]))


wms_DF <- left_join(saDF , tDF, by = c("x", "y"))
wms_XY <- wms_DF %>% select(x, y)
wms_DF %<>% select(-x,-y)

# POC flux:
pDF <- as_tibble(rasterToPoints(cor_Rasters[["poc"]]))

poc_DF <- pDF %>%  select(-x,-y)
poc_XY <- pDF %>%  select(x, y)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                                 Clustering                                               ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

gc()
# set a number of iteration for the clusters McQuaid et al. use 40
nb.iterations <- c(40,40,40)
# number of centers
nb.centers <- c(3,12,5)


group_name <- c("topo","wms","poc")

# nb.of clusters

# raw data
DFs <- list(topo_DF, wms_DF, poc_DF)
# cell positions
XYs <- list(topo_XY, wms_XY, poc_XY)
# rasters that will server as canvas
GridRasters <- list(cor_Rasters[["bbpi"]],cor_Rasters[["temperature"]], cor_Rasters[["poc"]])   

# for each group of layers calculate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# FOR EACH GROUPs OF VARIABLES
for (I in seq(group_name)) {
  
  # print the info
  print( 
  paste0(group_name[I]," Layers - ",nb.centers[I] ," centers", nb.iterations[I] ," iterations - - - - - - - - - - - - - - - - - - - - - -")
  )
  
  DF <- DFs[[I]]
  DFxy <- XYs[[I]]
  gridRaster <- GridRasters[[I]]
  
  # Calculate Average Silhouette Width:
  
  clara_clusters <-
    pamk(
      DF,
      krange = 2:nb.iterations,
      criterion = "asw",
      usepam = FALSE,
      scaling = FALSE,
      diss = inherits(DF, "Euclidean"),
      critout = TRUE
    )
  
  asw  <- as.data.frame(clara_clusters$crit)
  noclusters <- as.data.frame(as.numeric(c(1:nb.iterations)))
  results_asw <- cbind(noclusters, asw)
  names(results_asw) <- c("clusters", "asw")
  
  # Calculate Calinski-Harabasz index:
  clara_clusters <-
    pamk(
      DF,
      krange = 2:nb.iterations,
      criterion = "ch",
      usepam = FALSE,
      scaling = FALSE,
      diss = inherits(DF, "Euclidean"),
      critout = TRUE
    )
  
  ch <- as.data.frame(clara_clusters$crit)
  noclusters <- as.data.frame(as.numeric(c(1:nb.iterations)))
  results_ch <- cbind(noclusters, ch)
  names(results_ch) <- c("clusters", "ch")
  
  # Cluster based on ASW & CH:
  
  # Cluster based on ASW & CH:
  clara_clusters <- pamk(
    DF,
    krange = nb.centers[I], # number of centers adjusted manually
    criterion = "asw",
    usepam = FALSE,
    scaling = FALSE,
    diss = inherits(DF, "Euclidean"),
    critout = TRUE
  )
  
  clara_df <-
    data.frame(DFxy$x, DFxy$y, clara_clusters$pamobject$clustering)
  
  names(clara_df) <- c("x", "y", "cluster")
  coordinates(clara_df) <- ~ x + y
  clara_Ras <- rasterize(clara_df,  gridRaster)
  plot(clara_Ras[["cluster"]])
  
  # Export = = = = = = = = =
  
  write.csv(results_asw,
            file = paste0(NHADir, "/",group_name[I],"_asw.csv"),
            row.names = FALSE)
  write.csv(results_ch,
            file = paste0(NHADir, "/", group_name[I],"_ch.csv"),
            row.names = FALSE)
  writeRaster(
    clara_Ras,
    filename = paste0(NHADir, "/",  "clara_",group_name[I],".tif") ,
    format = "GTiff",
    progress = "text",
    overwrite = TRUE
  )
  
  
}# next group

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                          __Confusion Indices:                                            ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 

group_name <- c("topo","wms","poc")

varsList  <- list(c("bbpi","fbpi", "slope"),
                  c("temperature", "salinity"),
                  c("poc"))

nb.centers <- c(3,12,5)
nb.iters <- c(40,40,40)


# FOR EACH GROUPs OF VARIABLES

for (I in seq(group_name) ) {
  # to do 1 less
  
  print(varsList[I])
  
  var.I <-  raster::subset(var, varsList[[I]])
  
  #Convert rasters to vectors:
  
  tibble("ID" = 1:ncell(var.I)) -> data
  
  if( length(varsList[[I]]) > 1 ){
    
    vectors <- list()
    for (v in var.I@layers) {
      v %>%
        as.matrix() %>%
        as.vector() %>%
        bind_cols(data, .) -> data }
    # if there are one raster
  }else if(  length(varsList[[I]] ) == 1 ){ # when there is only 1 raster the loop doenst work
    
    var.I %>%
      as.matrix() %>%
      as.vector() %>%
      bind_cols(data, .) -> data
  }
  
  # Select variables, remove NAs and scale:
  
  data[data == 0] <- NA
  data_na <- na.omit(data)
  # remove the ID column
  data_norm <-
    (scale(select(data_na,-ID), scale = TRUE)) 
  data_norm <- as.data.frame(data_norm)
  
  # Run clustering: ----------------------------------------------------------------------------
  
  #
  fit <- kmeans(data_norm, nb.centers[I], iter.max = nb.iters[I])
  # the second argument is the number of clusters (here using the number determined through previous clustering)
  
  
  #Create matrix of input data & clustering output:
  
  
  m <- matrix(NA, nrow(data_norm), ncol = nrow(fit$centers))
  
  # calculate distance between each points and each centers
  # using Dist() function from R fast as it is a little faster than base R
  for (i in 1:nrow(data_norm)) {
    m[i, ] <- as.matrix(Dist(rbind(data_norm[i, ],  fit$centers)))[-1, 1]
  }
  
  mu <- matrix(NA, nrow(data_norm), ncol = nrow(fit$centers))
  
  #Calculate confusion index:
  for (i in 1:nrow(data_norm)) {
    mu[i, ] <- 1 / m[i, ] ^ 2 * 1 / sum(1 / m[i, ] ^ 2)
  }
  
  CI <- numeric(nrow(data_norm))
  
  for (i in 1:nrow(data_norm)) {
    CI[i] <-
      mu[i, order(mu[i, ], decreasing = TRUE)[2]] / mu[i, order(mu[i, ], decreasing =
                                                                  TRUE)[1]]
  }
  
  data_norm2 <- data.frame(data_norm, CI)
  matrix1 <- matrix(nrow = nrow(data), ncol = 2)
  temp1 <- cbind(data_na, data_norm2$CI)
  
  
  #Create table with correct dimensions:
  temp1 %>%  as.data.frame() %>%  as_tibble() %>%
    select(ID, `data_norm2$CI`) -> d
  
  matrix1 %>% as_tibble() %>%
    mutate(ID = 1:nrow(.)) %>%
    left_join(d, by = "ID") %>%
    rename(CI = `data_norm2$CI`) -> dd
  
  
  # a column for confidence indices
  temp1 %>%  as.data.frame() %>%  as_tibble() %>%
    select(ID, contains("$CI")) -> d
  
  # remake the raster and include the NA cells that had been removed
  matrix1 <- matrix(nrow = nrow(data), ncol = 2)
  # join the CI values with the cells that have one
  matrix1 %>% as_tibble() %>%
    mutate(ID = 1:nrow(.)) %>%
    left_join(d, by = "ID") %>%
    rename(CI = contains("$CI")) -> dd
  
  # make raster at the dimension of the input raster
  matrix3 <- matrix(dd$CI, nrow = nrow(var.I), ncol = ncol(var.I))
  
  #Plot and save output as raster:
  x <- raster(
    matrix3,
    xmn = -180,
    xmx = 180,
    ymn = -90,
    ymx = 90
  )
  plot(x, main = paste0(group_name[I], " Final classification"))
  plot(continents, add = T,col = "grey")
  
  crs(x) <- "+init=epsg:4326" # change coordinate system to WGS84
  writeRaster(x,
              filename = paste0(outputsDir, "/", "CI_Non-hier_", group_name[I], ".tif"))
  
  # remove results to ensure they are not reused by next iteration
  rm(data, x)
  
  
}













