title: "Regional benthic habitat classification of the South Atlantic"
author: "Tiago Gandra, Kirsty McQuaid, Oli Hogg & Nils Piechaud"
date: "April 2022"


# =================================================================================================================
##                                    PCAs 
# =================================================================================================================


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
# Add your own rasters in there. This shoud contains the layers (rasters and shapfile) downloaded in this repository: 

paste0(wd,"/outputs")  -> outputsDir
# in case you dont have it yet, create the directory
outputsDir %>%  dir.create()

paste0(outputsDir,"/NA") -> NADir
# in case you dont have it yet, create the directory
NADir %>%  dir.create()

paste0(outputsDir,"/HA") -> HADir
# in case you dont have it yet, create the directory
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


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                          __Principal Component Analysis__                                ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#PCA using "prcomp": - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

gc()
var.df <- as.data.frame(var)
var.df <- var.df[complete.cases(var.df), ]
head(var.df)
res.pca <- prcomp(var.df, scale = TRUE)
summary(res.pca)

#Set theme:
theme_set(
  theme_bw(base_size = 14) +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
)

# Screeplot:
plot(res.pca)
fviz_screeplot(res.pca, choice = 'eigenvalue', geom = 'line') +
  ggtitle(NULL) +
  geom_hline(yintercept = 1, linetype = 'dashed')

# export plot to output directory
ggsave(
  paste0(outputsDir, "/", "Level1_screeplot.jpg"),
  width = 13,
  heigh = 6,
  units = 'in',
  dpi = 150
)

# Loading plot:
fviz_pca_var(
  res.pca,
  axes = c(1, 2),
  # change 'c(,)' to desired PCs
  col.var = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE,
  title = ''
)

# export plot to output directory
ggsave(
  paste0(outputsDir, "/", "Level1_loadingplot.jpg"),
  dpi = 150,
  units = 'in',
  width = 7,
  height = 7
)

#PCA using "raster": - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

pca <- rasterPCA(var, spca = TRUE)
summary(pca$model)
plot(pca$map)

# Plot using plot() function:
pc <-
  stack(subset(pca$map, c(1:4))) # change 'c(1:4)' to desired PCs
plot(pc)

# writeRaster(pca$map, "Level1_pca.map.tif",   overwrite=TRUE) 

paste0(NADir,"/","pcRastersDir") -> pcrastersDir
# loop to export the rasters:
for (i in 1:nlayers(pca$map)) {
  pca$map[[i]] -> r
  writeRaster(r, paste0(pcrastersDir, "/", "Level1_", i, "_pca.map.tif"), overwrite = T)
}



# loop to export the rasters:
list.dirs(pcrastersDir,full.names = T) -> pcrasters
pc <- stack()
for (i in 1:nlayers(pca$map)) {
  
  raster( pcrasters[i] ) -> r
  
  addLayer(pc, r) -> pc
}

# Scaled maps:
pc_scale <-
  (pc - cellStats(pc, "min")) / (cellStats(pc, "max") - cellStats(pc, "min"))
pc_wgs <- projectRaster(pc_scale, crs = "+init=epsg:4326")
plot(pc_wgs)

# Plot using ggplot():
pc.p <- rasterToPoints(pc_scale)
pc.df <- data.frame(pc.p)

minlon <- min(pc.df$x - 0.2)
maxlon <- max(pc.df$x + 0.2)
minlat <- min(pc.df$y - 0.2)
maxlat <- max(pc.df$y + 0.2)

# convert to tall format
d <- reshape2::melt(pc.df, id = c('x', 'y'))

# save the table 

#
mapWorld <-
  borders("world",
          colour = "gray50",
          fill = "gray50",
          alpha = .7)

ggplot() + geom_raster(data = d, aes(y = y, x = x, fill = value)) + theme_bw() +
  facet_wrap(~ variable, ncol = 2) +
  coord_equal(xlim = c(minlon, maxlon),
              ylim = c(minlat, maxlat)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_distiller(palette = "Spectral") +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 8),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    line = element_blank()
  ) +
  mapWorld

# export plot to output directory
ggsave(
  paste0(outputsDir, "/", 'Level1_scaled map_ggplot.jpg'),
  dpi = 150,
  units = 'in',
  width = 13,
  height = 7.5
)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                             #Calculate correlation between environmental variables and PCs:              ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 


#Calculate correlation between environmental variables and PCs:

var_pca <- stack(var, pc)
pca.df <- as.data.frame(var_pca)
pca.df <- pca.df[complete.cases(pca.df), ]
c <- cor(pca.df, method = "pearson")
write.table(c, file = 'Level1_correlations.csv', dec = '.', sep = ';')

# Plot correlogram:
png(
  height = 500,
  width = 700,
  file = paste0(outputsDir, "/", "Level1_corplot.png"),
  type = "cairo"
)
corrplot(c, method = 'square', type = 'lower')
dev.off()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                            #__Clustering Principal Components__                                          ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#Calculate Calinski-Harabasz index:  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# turn the principla components rasters into a table
df <- as.data.frame(scale(pc))
df <- as.data.frame(df[complete.cases(df), ])

# Change 2:40 to desired number of iterations:
ci <- data.frame(1, NA)
names(ci) <- c('i', 'cii')

# try up to how many numbers of clusterscc
nbiterations <- 40


for (i in c(2:nbiterations)) {
  # print how much there is left
  paste0(((i / nbiterations) * 100) %>% round() , "% done") %>% print()
  # run clustering
  c <- kmeans(df, i,  iter.max = 10, nstart = 1)
  cii <- calinhara(df, c$cluster, cn = i)
  ci <- rbind(ci, data.frame(i, cii))
}

ci %>% arrange(-cii)


#Plot results, changing intercept results to peaks in CH:
ggplot(ci, aes(x = i, y = cii / 10000)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.3) +
  theme_bw() + xlab('Clusters') + ylab('c-Harabasz Index') +
  geom_vline(xintercept = 5, linetype = 2) +
  geom_vline(xintercept = 4, linetype = 2) +
  geom_vline(xintercept = 6, linetype = 2)

ggsave(
  paste0(outputsDir, "/", "Level1_CHindex.jpg"),
  width = 13 ,
  height = 6,
  units = 'in',
  dpi = 200
)

#Calculate Average Silhouette Width:  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

pc.df <- pc.df[complete.cases(pc.df), ]
dim(pc.df)
head(pc.df)

pc.pamk <-
  pc.df[, c(3:6)] # Change 'c(3:6)' to capture columns with PC data
clara <- pamk(
  pc.pamk,
  krange = 2:40,
  criterion = "asw",
  usepam = FALSE,
  scaling = FALSE,
  diss = inherits(pc.pamk, "Euclidean"),
  critout = TRUE
)

#Cluster analysis:

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
clara.rast <- rasterFromXYZ(clara.out)
crs(clara.rast) <- "+init=epsg:4326"
plot(clara.rast)
writeRaster(clara.rast,  paste0(outputsDir, "/", "Level1.tif"))

#__Boxplots of cluster outputs__

var.crop <- crop(var, clara.rast)

s <- stack(var.crop, clara.rast)

s.df <- as.data.frame(s)
s.df <- s.df[complete.cases(s.df), ]

unique(s.df$zone)
s.m <- melt(s.df, id = 'zone')

#Define number of colors for plot:
nb.cols <- max(s.df$zone, na.rm = T)
mycolors <- colorRampPalette(brewer.pal(4, "YlGnBu"))(nb.cols)

# take a smaller subset of the dataset to do the plot - Change fraction to make the plot more faithful to the data
# !! This step can results in crach as it consumes a lot of memory !!

s.m %>% sample_frac(0.01) -> s.m2

#Plot boxplots:
ggplot(s.m2, aes(
  x = factor(zone),
  y = value,
  fill = factor(zone)
)) +
  geom_boxplot(outlier.size = .1, lwd = .1) +
  facet_wrap( ~ variable, scales = 'free', ncol = 4) +
  theme_bw() +
  scale_fill_manual(values = mycolors) +
  theme(
    legend.position = "none",
    text = element_text(family = "serif"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  xlab("Habitat class") + ylab(NULL)

ggsave(
  paste0(outputsDir, "/", 'Level1_boxplots.jpg'),
  width = 16,
  height = 6.5,
  unit = 'in',
  dpi = 300
)
gc()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                        #__Decision tree__                                                ========
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#__Decision tree__ 

s<-stack(var.crop,clara.rast)
s.df<-as.data.frame(s)
s.df<-s.df[complete.cases(s.df),]
levels(s.df$zone) <- c("Class 1", "Class 2", "Class 3", "Class 4")
s.df$zone<-factor(s.df$zone) 

#Decision tree:

bfit <- rpart(zone~., method="class", data=s.df) 
printcp(bfit)  
plotcp(bfit) 
summary(bfit) 

# Save the plots

jpeg(paste0(outputsDir,"/",'Level1_decisiontree.jpg'), width=13, height = 7, units = 'in', res = 300)
rpart.plot(bfit, type=0, extra = 100, under = TRUE, nn=FALSE, family = 'serif', box.palette = 'Grays')
dev.off()
