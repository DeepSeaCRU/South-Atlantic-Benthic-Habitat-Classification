title: "Regional benthic habitat classification of the South Atlantic"
author: "Tiago Gandra, Kirsty McQuaid, Oli Hogg & Nils Piechaud"
date: "April 2022"

  
# =================================================================================================================
##                                    Level 2 clustering 
# =================================================================================================================

#Level 2 of the hierarchical classification is nested within Level 1. Level 2 clusters the environmental data for each output habitat class of Level 1 in isolation. Code for plots and figures to support these analyses can be found in Level 1 script.



# list packages that will be needed
pacakges <- c("sp","raster","rworldmap", "rgdal","factoextra","RStoolbox", "corrplot", "fpc","RColorBrewer", "rpart", "rpart.plot" ,"NbClust",  "psych"  , "reshape2", "tidyverse", "magrittr")

# install those packages that you dont have yet
install.packages( 
  setdiff(pacakges,installed.packages()[,1] )
)


#Load packages:

library(raster)
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
    "sal",
    "tmp",
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

 
#Load output of Level 1 clustering:  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
paste0(outputsDir,"/","Level1.tif") %>% raster() -> zones
names(zones) <- 'zone' # output of level 1 clustering
 

#Crop environmental variables to extent of each habitat class in Level 1:
 
# Crop raster stack to extent of clara output:
var %<>%  crop(zones)

# plot each zone
zones.values <- values(zones)  %>%  na.omit()%>%  unique() 

ZonesRasters <- as.list(zones.values)
names(ZonesRasters) <- zones.values

 for (v in zones.values) {
  zone.v <-  zones == v
  
  zone.v[zone.v < 1] <- NA # this new raster has values 0-1 (i.e. binomial: zone 2 yes or no), so select only those values = 1
  # Crop raster stack to extent of each habitat class:
  
  zone.v.var <- mask(var , zone.v)
  
  ZonesRasters[[v]] <- zone.v.var
  
  # plot the extent of each zone -- just add # to cancel rows
      #plot(zone.v, main = paste0("Zone ",v))
       plot(zone.v.var$sal, main = paste0("Zone ",v, " - Salinity"))
       plot(continents, add = T, col = "grey")
}

names(ZonesRasters) <- paste0("zone",zones.values)




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
##                                    __Clustering of Class 1 to 4                                            =====
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 


names(ZonesRasters) -> ZONES

# set the number of iterations in the kmean clustering
nb.iterations <- 40

# how many PCs ?
tibble(ZONES, nb.pcs = c(3,4,4,4)) -> paramters

# Average Silhouette Width:
  paramters %<>% 
    mutate(sil.wdth.left = c(3,3,3,3), 
           sil.wdth.right = c(5,6,6,6))




for (I in seq(ZONES)) {
  
  # select parameters 
  paramters %>% filter(ZONES == ZONES[I]) -> p
  
  
  
  gc()  
  #Principal Component Analysis:
  pca <- rasterPCA( ZonesRasters[[ZONES[I]]] ,spca=TRUE)
  summary(pca$model)
  
  #Clustering Principal Components:
  
  # Calculate correlation between environmental variables and PCs:
  pc <- stack(subset(pca$map,c( 1:p$nb.pcs ))) # change 1:3 based on desired number of PCs
  zoneIvar_pca <- stack(ZonesRasters[[ZONES[I]]], pc)
  pca.df <- as.data.frame(zoneIvar_pca)
  pca.df <- pca.df[complete.cases(pca.df),]
  ccor <- cor(pca.df, method="pearson")
  
  write.table(
    ccor,
    file = paste0(HADir, "/", 'Level2_correlations_zone', I,  '.csv'),
    dec = '.',
    sep = ';'
  )
   
  
  # Calculate Calinski-Harabasz index:
  pc_scale <- (pc-cellStats(pc,"min"))/(cellStats(pc,"max")-cellStats(pc,"min"))
  df <- as.data.frame(scale(pc))
  df<-as.data.frame(df[complete.cases(df),])
  
  # table to store results of the clustering
  ci <- data.frame(1,NA)
  names(ci) <- c('i','cii')
  
  # Kmeans
  for(i in c(2:nb.iterations)){
    c<-kmeans(df, i)
    cii<-calinhara(df,c$cluster,cn=i)
    ci<-rbind(ci,data.frame(i,cii))
  }
  
  
  # Calculate Average Silhouette Width:
  pc.p <- rasterToPoints(pc_scale)
  pc.df <- data.frame(pc.p)
  head(pc.df)
  pc.pamk <- pc.df %>% select(-x, - y) # only keep the PCs columns
  clara <-
    pamk(
      pc.pamk,
      krange = 2:nb.iterations,
      criterion = "asw",
      usepam = FALSE,
      scaling = FALSE,
      diss = inherits(pc.pamk, "Euclidean"),
      critout = TRUE
    )
  
  
  # Cluster analysis:
  clara <- pamk(
    pc.pamk,
    krange = 4,
    criterion = "asw",
    usepam = FALSE,
    scaling = FALSE,
    diss = inherits(pc.pamk, "Euclidean"),
    critout = TRUE
  ) # change 4 above to desired number of clusters
  
  
  
  
  clara.out <-
    data.frame(pc.df$x, pc.df$y, clara$pamobject$clustering)
  zone <- rasterFromXYZ(clara.out)
  crs(zone) <- "+init=epsg:4326"
  plot(zone, main = paste0("cluster Analysis - ",ZONES[I]))
  writeRaster(zone, paste0(HADir, "/", "Level2_Class", I, ".tif") , overwrite = TRUE)
  
  # delete output R object
  rm(pca,ccor,zone, c)
  
  
}
  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
##                                    __Crop all Level 2 outputs to South Atlantic__                           ====
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# import mask of South Atlantic  
SAmask <- readOGR( paste0(rasterDir , "/", "Mask_SouthAtlantic.shp") )
  
  Level2Rasters <-
    HADir %>% list.files(pattern = "Level2_Class", full.names = T)
  Level2Rasters_names <-
    HADir %>% list.files(pattern = "Level2_Class", full.names = F) %>%
    str_remove(pattern = ".tif") 


for (i in seq(Level2Rasters)) {
  
  
  print(Level2Rasters_names[i])
  
  Class <- raster(Level2Rasters[i])
  # crop with South Atlantic Mask
  Class.crop <- crop(Class, SAmask)
  
  # 
  Class <- mask(Class.crop, SAmask)
  
  # export croped raster
  writeRaster(Class,  paste0(HADir,"/",Level2Rasters_names[i],"_SAtlant.tif" ))
  
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
##                                    Level 3 clustering                                                       ====
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  
  
  
## 
  
  #Level 3 of the hierarchical classification is nested within Level 2.
  # Level 3 clusters the environmental data for each output habitat class of Level 2 in isolation. 
  #Code to support Level 3 clustering follows the same pattern as code for Level 2, 
  #the only difference being the use of Level 2 cluster outputs to crop environmental data (instead of Level 1). 
  # For plots and figures to support these analyses, see Level 1 script.  
  
  

 

