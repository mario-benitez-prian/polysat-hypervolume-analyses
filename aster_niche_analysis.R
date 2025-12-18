library(readxl)
library(tidyverse)#data management
library(raster)#to read climate data
library(dplyr)
library(hypervolume)##For hypervolume analysis
library(ggplot2)

getwd()
setwd("C:/Users/mario/Desktop/Supplemmentary_Material_Aster_amellus/R_project")

### Reading the location data

amg <- read_excel("../scripts_data/aster_contact_zone.xlsx")
summary(amg)
str(amg)
amg

### Reading and extracting from raster

path.wc = "Climate data"
s.wc = list.files(path=paste(path.wc), pattern='tif', full.names=TRUE)
str(s.wc)
s.wc = stack(s.wc)

#### Extracting worldclim values

g.clm.wc=extract(s.wc, cbind(amg$Lon, amg$Lat), cellnumbers=T)
head(g.clm.wc)

g.clm.wc= data.frame(g.clm.wc)

colnames(g.clm.wc)

newcol.names <- c("bio1","bio10","bio11","bio12","bio13","bio14","bio15",
                  "bio16","bio17","bio18","bio19",
                  "bio2","bio3","bio4","bio5", 
                  "bio6","bio7","bio8","bio9")

colnames(g.clm.wc)[2:20] <- newcol.names

colnames(g.clm.wc)

### Merge locations with climatic variables 

sh.wc= cbind(amg,g.clm.wc) 
head(sh.wc)


### PCA 

pca = prcomp((sh.wc[10:28]), scale = TRUE)
pca
summary(pca)
pca$rotation

### Selecting two first axes 

pca.axis <- data.frame(PC1= pca$x[,1], PC2=pca$x[,2])

amg.pca <- data.frame(cbind(sh.wc, pca.axis))###combining the principal components in the dataframe

head(amg.pca)


### PCA scatterplots

biplot(pca,
       cex = 1)


### Diploid hypervolume 

library(alphahull)

dip = amg.pca %>% filter(Ploidy == "2x") 

dip.hv = hypervolume(data = cbind(dip$PC1, dip$PC2),
                    method = "svm", name = "diploid svm")
dip.hv

plot(dip.hv, 
     col = c(rgb(red = 1, green = 0, blue = 0, alpha = 0.4), rgb(red = 0, green = 0, blue = 1, alpha = 0.4)),
     show.random=FALSE, show.density=TRUE,show.data=TRUE,
     show.contour=TRUE,
     contour.lwd=1, 
     contour.type='raster', 
     contour.raster.resolution=32,
     show.centroid=TRUE, cex.centroid=FALSE
)



### Hexaploid hypervolume

hex = amg.pca %>% filter(Ploidy == "6x") 

hex.hv = hypervolume(data = cbind(hex$PC1, hex$PC2),
                     method = "svm", name = "hexaploid svm")
hex.hv

plot(hex.hv,
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.4),
     show.random=FALSE, show.density=TRUE,show.data=TRUE,
     show.contour=TRUE,
     contour.lwd=1, 
     contour.type='raster', 
     contour.raster.resolution=32,
     show.centroid=TRUE, cex.centroid=FALSE)

###Hypervolume set - to get union, intersection and unique components from volumes

(h.set = hypervolume_set(dip.hv, hex.hv, check.memory = F))

### Overlap statistics - Jaccard, Sorensen, unique components

hypervolume_overlap_statistics(h.set)


### Get p-values for overlap statistics

qsample1 = dip[sample(1:nrow(dip), 15),]
qsample2 = hex[sample(1:nrow(hex), 16),]


hv_combine = hypervolume(rbind(cbind(qsample1$PC1, qsample1$PC2), cbind(qsample2$PC1, qsample2$PC2)))
bootstrap_path2 = hypervolume_resample("hypervolume_bootstrap2",
                                      hv_combine,
                                      method = "bootstrap",
                                      n = 100,
                                      points_per_resample = 30,
                                      cores = 1)

overlap_test <- hypervolume_overlap_test(dip.hv, hex.hv, bootstrap_path2, alternative = "one-sided", bins = 100, cores = 1)

### get centroid distance 

hypervolume_distance(dip.hv, hex.hv, type = "centroid",
                     num.points.max = 1000, check.memory = TRUE)

### Figure 6

hv_list = hypervolume_join(dip.hv,hex.hv)

plot(hv_list, 
     #colors = c("#e0473f", "#729be8"),
     col = c(rgb(red = 1, green = 0, blue = 0, alpha = 0.3), rgb(red = 0, green = 0, blue = 1, alpha = 0.3)),
     show.random=FALSE, show.density=TRUE,show.data=TRUE,
     show.contour=TRUE,
     contour.lwd=2, 
     contour.type='raster', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=2, 
     contour.kde.level=1e-04,
     contour.raster.resolution=32,
     show.centroid=TRUE, cex.centroid=4
)
