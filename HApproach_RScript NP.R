---
title: "Regional benthic habitat classification of the South Atlantic: Hierarchical approach"
author: "Tiago Gandra, Kirsty McQuaid, Oli Hogg & Nils Piechaud"
date: "April 2022"
output: html_document
---

  
  
 
  
  
  
#Load packages:
#```{r}  
library(raster)
library(factoextra)
library(RStoolbox)
#library(reshape2)
library(corrplot)
library(fpc)
library(RColorBrewer)
# library(dplyr)
library(rpart)
library(rpart.plot)
#library(dplyr)
library(NbClust)
library(sp)
library(psych)
library(tidyverse)
library(rgdal)
library(magrittr)
#library(ggtext)  
#```

## Level 1 clustering

#__Load environmental variables__

#Set working directory to where environmental input data are saved:
#```{r}
setwd("C:/Users/kmcquaid/Desktop/Nils/inputs")
#```

#Load data:
#```{r}
fbpi <- raster("GEBCO2020_FBPI_10km_WGSfromMOL.tif")
bbpi <- raster("GEBCO2020_BBPI_10km_WGSfromMOL.tif")
slope <- raster("GEBCO2020_Slope_10km_WGSfromMOL.tif")
poc <- raster("2006-2015_fluxXYZ112deg.tif")
currvel <- raster("Present.Benthic.Mean.Depth.Current.Velocity.Mean.tif.BOv2_1.tif")
dissox <- raster("Present.Benthic.Mean.Depth.Dissolved.oxygen.Mean.tif")
nit <- raster("Present.Benthic.Mean.Depth.Nitrate.Mean.tif")
phos <- raster("Present.Benthic.Mean.Depth.Phosphate.Mean.tif")
sal <- raster("BioO_sal_10km_WGS84.tif")
sil <- raster("Present.Benthic.Mean.Depth.Silicate.Mean.tif")
tmp <- raster("BioO_temp_10km_WGS84.tif")
depth <- raster("GEBCO2020_depth_10km_WGSfromMOL.tif")
#```

#Set working directory to results folder:
#```{r}
setwd("../outputs/HA")
#```

#Make raster stack:
#```{r}
# Crop all variables to same extent
stack1 <- stack(poc, currvel, dissox, nit, phos, sal, sil, tmp)
stack2 <- stack(fbpi, bbpi, slope)
stack2.crop <- crop(stack2, stack1)
# Stack all variables
var <- stack(stack1, stack2.crop)
names(var)
names(var)<-c('poc','currvel','dissox','nitrate','phosphate','salinity','silicate','temp','fbpi', 'bbpi', 'slope')
#```

#Create mask:
#```{r}
mask_<-depth<=0
mask <- crop(mask_, var )
plot(mask)
#```

#__Principal Component Analysis__
                             
#PCA using "prcomp": # ======================================================================================
#```{r}
gc()
var.df<-as.data.frame(var)
var.df<-var.df[complete.cases(var.df),]
head(var.df)
res.pca<-prcomp(var.df,scale = TRUE)
summary(res.pca)
#```

#Plot PCA results:
#```{r}
#Set theme:
theme_set(
  theme_bw(base_size = 14)+
    theme(text=element_text(family="Arial"),
          plot.title = element_text(hjust = 0.5, face='bold',size=14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

# Screeplot:
plot(res.pca)
fviz_screeplot(res.pca, choice='eigenvalue', geom='line')+
  ggtitle(NULL)+
  geom_hline(yintercept=1, linetype='dashed')
ggsave('Level1_screeplot.jpg', width=13, heigh=6, units='in',dpi=150)

# Loading plot:
fviz_pca_var(res.pca, axes = c(1,2),    # change 'c(,)' to desired PCs
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     
             title='')
ggsave("Level1_loadingplot.jpg", dpi=150, units='in', width=7, height=7)
#```

#PCA using "raster": # ======================================================================================
#```{r}
pca<-rasterPCA(var,spca=TRUE)
summary(pca$model)
plot(pca$map)
#```

#Plot PCA maps:
#```{r}
# Plot using plot() function:
pc<-stack(subset(pca$map,c(1:4))) # change 'c(1:4)' to desired PCs
plot(pc)
writeRaster(pca$map, "Level1_pca.map.tif")

# Scaled maps:
pc_scale<-(pc-cellStats(pc,"min"))/(cellStats(pc,"max")-cellStats(pc,"min"))
pc_wgs<-projectRaster(pc_scale,crs="+init=epsg:4326")
plot(pc_wgs)

# Plot using ggplot():
pc.p <- rasterToPoints(pc_scale)
pc.df <- data.frame(pc.p)

minlon<-min(pc.df$x-0.2)
maxlon<-max(pc.df$x+0.2)
minlat<-min(pc.df$y-0.2)
maxlat<-max(pc.df$y+0.2)

d<-melt(pc.df,id=c('x','y'))
mapWorld <- borders("world", colour="gray50", fill="gray50", alpha=.7)

ggplot()+geom_raster(data=d,aes(y=y,x=x, fill=value))+theme_bw()+
  facet_wrap(~variable, ncol=2)+
  coord_equal(xlim = c(minlon, maxlon),ylim = c(minlat, maxlat))+
  xlab(NULL)+ylab(NULL)+
  scale_fill_distiller(palette = "Spectral")+
  theme(legend.position="none",
        strip.text.x = element_text(size = 8),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        line=element_blank())+
mapWorld

ggsave('Level1_scaled map_ggplot.jpg', dpi=150, units='in', width=13, height=7.5)
#```

#Calculate correlation between environmental variables and PCs: # ======================================================================================
#```{r}
var_pca<-stack(var, pc)
pca.df<-as.data.frame(var_pca)
pca.df<-pca.df[complete.cases(pca.df),]
c<-cor(pca.df, method="pearson")
write.table(c, file='Level1_correlations.csv', dec = '.', sep=';')

# Plot correlogram:
png(height=500, width=700, file="Level1_corplot.png", type = "cairo")
corrplot(c,method='square',type='lower')
dev.off()
#```

#__Clustering Principal Components__ # ======================================================================================


#Calculate Calinski-Harabasz index:
#```{r}
df <- as.data.frame(scale(pc))
df<-as.data.frame(df[complete.cases(df),])

# Change 2:40 to desired number of iterations:
ci=data.frame(1,NA)
names(ci)=c('i','cii')
nbiterations <- 40


for(i in c(2:nbiterations)){
  # print how much there is left
  paste0( ((i/nbiterations)*100) %>% round() ,"% done" )
  
  c<-kmeans(df, i)
  cii<-calinhara(df,c$cluster,cn=i)
  ci<-rbind(ci,data.frame(i,cii))
}

ci%>%arrange(-cii)

#Plot results, changing intercept results to peaks in CH:
ggplot(ci,aes(x=i,y=cii/10000))+geom_point(size=.5)+geom_line(size=.3)+
  theme_bw()+xlab('Clusters')+ylab('Calinski-Harabasz Index')+
  geom_vline(xintercept = 5, linetype=2)+
  geom_vline(xintercept = 4, linetype=2)+
  geom_vline(xintercept = 6, linetype=2)

ggsave("Level1_CHindex.jpg", width = 13 , height = 6, units = 'in', dpi = 200)
#```

#Calculate Average Silhouette Width: # ======================================================================================

#```{r}
pc.df<-pc.df[complete.cases(pc.df),]
dim(pc.df)
head(pc.df)

pc.pamk <- pc.df[,c(3:6)] # Change 'c(3:6)' to capture columns with PC data
clara <- pamk(pc.pamk, krange=2:40,criterion="asw", usepam=FALSE,
              scaling=FALSE, diss=inherits(pc.pamk, "Euclidean"),
              critout=TRUE)
#```

#Cluster analysis:
#```{r}
clara <- pamk(pc.pamk, krange=4,criterion="asw", usepam=FALSE,
              scaling=FALSE, diss=inherits(pc.pamk, "Euclidean"),
              critout=TRUE) # change 4 above to desired number of clusters

clara.out <- data.frame(pc.df$x, pc.df$y, clara$pamobject$clustering)
clara.rast <- rasterFromXYZ(clara.out)
crs(clara.rast) <- "+init=epsg:4326"
plot(clara.rast)
writeRaster(clara.rast, "Level1.tif")
#```

#__Boxplots of cluster outputs__ # ======================================================================================


#Prepare environmental variables:
#```{r}
var.crop <- crop(var, clara.rast) 
names(clara.rast)<-'zone' 
s <- stack(var.crop, clara.rast) 
names(s)
names(s)<-c('POC','Current velocity','Dissolved oxygen','Nitrate','Phosphate','Salinity','Silicate','Temperature','FBPI','BBPI', 'Slope', 'zone')

s.df<-as.data.frame(s)
s.df<-s.df[complete.cases(s.df),]

unique(s.df$zone)
s.m<-melt(s.df,id='zone')
#```

#Define number of colors for plot:
#```{r}
nb.cols <- max(s.df$zone,na.rm=T)
mycolors <- colorRampPalette(brewer.pal(4, "YlGnBu"))(nb.cols)
#```

#Plot boxplots:
#```{r}
ggplot(s.m,aes(x=factor(zone),y=value,fill=factor(zone)))+
  geom_boxplot(outlier.size=.1,lwd=.1)+
  facet_wrap(~variable,scales='free',ncol=4)+
  theme_bw()+
  scale_fill_manual(values=mycolors)+
  theme(
    legend.position="none",
    text = element_text(family = "serif"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))+
  xlab("Habitat class")+ylab(NULL)

ggsave('Level1_boxplots.jpg', width=16, height=6.5, unit='in', dpi=300)
#```

#__Decision tree__ # ======================================================================================


#Prepare environmental variables:
#```{r}
s<-stack(var.crop,clara.rast)
names(s)<-c('poc','currvel','dissox','nitrate','phosphate','salinity','silicate','temp','fbpi','bbpi', 'slope', 'zone')
s.df<-as.data.frame(s)
s.df<-s.df[complete.cases(s.df),]
levels(s.df$zone) <- c("Class 1", "Class 2", "Class 3", "Class 4")
s.df$zone<-factor(s.df$zone)
#```

#Decision tree:
#```{r}
bfit <- rpart(zone~., method="class", data=s.df) 
printcp(bfit)  
plotcp(bfit) 
summary(bfit) 
#```

#Save plot:
#```{r}
jpeg('Level1_decisiontree.jpg', width=13, height = 7, units = 'in', res = 300)
rpart.plot(bfit, type=0, extra = 100, under = TRUE, nn=FALSE, family = 'serif', box.palette = 'Grays')
dev.off()
#```

#__Confusion Index__ # ======================================================================================


#calculating confusion index for Level 1 clustering of the hierarchical approach. Confusion index values range from 0 to 1, with values approaching 1 indicating greater uncertainty in the clara clustering solution between two or more clusters.

#Convert rasters to vectors:
#  ```{r}
unstack(var)
names(var)

raster1 <- var$fbpi
raster2 <- var$bbpi
raster3 <- var$slope
raster4 <- var$poc
raster5 <- var$currvel
raster6 <- var$dissox
raster7 <- var$nitrate
raster8 <- var$phosphate
raster9 <- var$salinity
raster10 <- var$silicate
raster11 <- var$temp

matrix1 <- as.matrix(raster1)
matrix2 <- as.matrix(raster2)
matrix3 <- as.matrix(raster3)
matrix4 <- as.matrix(raster4)
matrix5 <- as.matrix(raster5)
matrix6 <- as.matrix(raster6)
matrix7 <- as.matrix(raster7)
matrix8 <- as.matrix(raster8)
matrix9 <- as.matrix(raster9)
matrix10 <- as.matrix(raster10)
matrix11 <- as.matrix(raster11)

vector1 <- as.vector(matrix1)
vector2 <- as.vector(matrix2)
vector3 <- as.vector(matrix3)
vector4 <- as.vector(matrix4)
vector5 <- as.vector(matrix5)
vector6 <- as.vector(matrix6)
vector7 <- as.vector(matrix7)
vector8 <- as.vector(matrix8)
vector9 <- as.vector(matrix9)
vector10 <- as.vector(matrix10)
vector11 <- as.vector(matrix11)
#  ```

#Create dataframe, remove NAs and scale:
#  ```{r}
data <- data.frame(vector1, vector2, vector3, vector4, vector5, vector6, vector7, vector8, vector9, vector10, vector11)
data[1:9331200, "ID"] <- c(1:9331200)
data[data==0] <- NA
data_na <- na.omit(data)
data_norm <- (scale(data_na[,-12], scale=TRUE)) ######## What is the number -4 for?? In Oli's code this value was his number of variables + 1. So for now have just made it my no. vars + 1
data_norm <- as.data.frame(data_norm)
#```

#Run Principal Component Analysis and clustering:
#  ```{r}
# PCA:
data_PC <- prcomp(data_norm, scale=FALSE)
eigen <- data_PC$sdev^2
eigen # Examine eigenvalues to identify number of PCs with eigenvalues > 1
dataPC_rotated <- psych::principal(data_norm, rotate="varimax", nfactors=4, scores=TRUE) # Change nfactors to desired number of PCs
data_PCs <- dataPC_rotated$scores

# Clustering:
fit <- kmeans(data_PCs, 4, iter.max=40)

#Create matrix of input data & clustering output:
#  ```{r}
m <- matrix(NA, nrow(data_PCs), ncol=nrow(fit$centers))

for(i in 1:nrow(data_PCs)){
  m[i,]=as.matrix(dist(rbind(data_PCs[i,],  fit$centers)))[-1,1]
}

mu <- matrix(NA, nrow(data_PCs), ncol=nrow(fit$centers))
#```

#Calculate confusion index:
#  ```{r}
for(i in 1:nrow(data_PCs)){
  mu[i,]=1/m[i,]^2*1/sum(1/m[i,]^2)
}

CI <- numeric(nrow(data_PCs))

for(i in 1:nrow(data_PCs)){
  CI[i] <- mu[i,order(mu[i,], decreasing=TRUE)[2]]/mu[i,order(mu[i,], decreasing=TRUE)[1]]
}

data_PCs2 <- data.frame(data_PCs, CI)
matrix1 <- matrix(nrow=nrow(data ), ncol=2)
temp1 <- cbind(data_na, data_PCs2$CI)
#  ```

#Create table with correct dimensions:
#  ```{r}
temp1 %>%  as.data.frame() %>%  as_tibble() %>%
 select(ID,`data_PCs2$CI`) -> d
  matrix1 %>% as_tibble() %>% mutate(ID = 1:nrow(.)) %>%  left_join(d) %>% rename(CI = `data_PCs2$CI`) -> dd
matrix3 <- matrix(dd$CI, nrow=2160, ncol=4320)
#```

#Plot and save output as raster:
#  ```{r}
Level1_CI <- raster(matrix3, xmn=-180, xmx=180, ymn=-90, ymx=90)
plot(Level1_CI)
crs(Level1_CI) <- "+init=epsg:4326"
writeRaster(Level1_CI, filename="Level1_CI.tif")

#__Crop outputs to South Atlantic study area:__

#Load mask of South Atlantic:
#  ```{r}
setwd("../../inputs")
SAmask <- readOGR("Mask_SouthAtlantic.shp")
#  ```

#Crop classification Level 1:
#  ```{r}
# level 1
level1 <- clara.rast
level1.crop <- crop(level1, SAmask)
level1 <- mask(level1.crop, SAmask)
plot(level1)
writeRaster(level1, "../outputs/HA/Level1_SAtlant.tif")
#  ```

#Crop confidence map Level 1:
#  ```{r}
Level1_CI.crop <- crop(Level1_CI, SAmask)
Level1_CI <- mask(Level1_CI.crop, SAmask)
plot(Level1_CI)
writeRaster(Level1_CI, "../outputs/HA/Level1_CI_SAtlant.tif")
#  ```