library(polysat)
library(ggplot2)
library(ggfortify)
library(hrbrthemes)

# DATA MATRIX SETTING 

# Set the scripts_data directory path
setwd("C:/Users/mario/Desktop/Aster amellus/Supplemmentary_Material/scripts_data")

data <- read.GeneMapper("data_polysat.txt")
summary(data)
Samples(data)
Loci(data)
find.missing.gen(data)

#Loading PopNames, PopInfo + Ploidies tables

popNamesTable <- read.table("pop_names.txt", header = TRUE)
popInfoTable <- read.table("PopInfo_Ploidies.txt", header= TRUE)

#Setting auxiliar info (population names, population IDs, population ploidies 
#and nucleotidic repetitions)

Description(data) <- "Aster amellus microsatellite data"
PopNames(data) <- popNamesTable[[1]]
PopInfo(data) <- popInfoTable[[2]]
Ploidies(data) <- popInfoTable[[3]]
Usatnts(data) <- c(2, 3, 4, 2, 3, 3, 2)

#Saving final matrix
save(data, file="data_aster.RData")

#                   SETTING DATA SETS FOR OTHER PROGRAMS 

write.Structure(data, ploidy = 6, file="data_structure.txt")
write.GenoDive(data, file="data_genoDive.txt", digits = 3)

#CREATING DATA SUBSETS FOR STRUCTURE ANALYSES 

#Creating subsets by ploidy
#You have to index ploidies by sample 
ploidia_sample <- reformatPloidies(data, output = "sample", na.rm = TRUE, 
                                   erase = FALSE)
Ploidies(ploidia_sample)
diploid_pops <- Samples(ploidia_sample, ploidies = 2)
hexaploid_pops <- Samples(ploidia_sample, ploidies = 6)

#Creating subsets by population
South_pops <- Samples(data, populations = c("Pop5", "Pop7", "Pop11", "Pop17", 
                                            "Pop18", "Pop20", "Pop25", "Pop26"))

NW_pops <- Samples(data, populations = c("Pop2", "Pop4", "Pop6", "Pop8", "Pop9",
                                         "Pop13", "Pop15", "Pop16", "Pop21",
                                         "Pop23", "Pop29", "Pop30", "Pop31",
                                         "Pop32", "Pop33"))

East_pops <- Samples(data, populations = c("Pop1", "Pop3", "Pop10", "Pop12", 
                                           "Pop19", "Pop22", "Pop24", "Pop27", 
                                           "Pop28"))

kostalov_pop <- Samples(data, populations = "Pop21")

mixed_pop <- Samples(data, populations = c("Pop7", "Pop20"))

#Get the number of samples by subset, you need this for structure analyses config
length(South_pops)
length(NW_pops)
length(East_pops)
length(kostalov_pop)
length(mixed_pop)
length(diploid_pops)
length(hexaploid_pops)

#Writing subsets to Structure data format
write.Structure(data, ploidy = 6, file = "NW_pops_structure.txt", samples = NW_pops, writepopinfo = TRUE, missingout = -9)
write.Structure(data, ploidy = 6, file = "South_pops_structure.txt", samples = South_pops, writepopinfo = TRUE, missingout = -9)
write.Structure(data, ploidy = 6, file = "East_pops_structure.txt", samples = East_pops, writepopinfo = TRUE, missingout = -9)
write.Structure(data, ploidy = 2, file = "diploid_pops_structure.txt", samples = diploid_pops, writepopinfo = TRUE, missingout = -9)
write.Structure(data, ploidy = 6, file = "hexaploid_pops_structure.txt", samples = hexaploid_pops, writepopinfo = TRUE, missingout = -9)
write.Structure(data, ploidy = 6, file = "kostalov_pop_structure.txt", samples = kostalov_pop, writepopinfo = TRUE, missingout = -9)
write.Structure(data, ploidy = 6, file = "mixed_pop_structure.txt", samples = mixed_pop, writepopinfo = TRUE, missingout = -9)


#Data to SPageDi format 
write.SPAGeDi(data, file="dataSpagedi.txt", digits = 3)



#                   ANALYSES 

#Create matrix with bruvo distances 
bruvo_dist <- meandistance.matrix(data, samples = Samples(data), 
                                  loci = Loci(data), all.distances = TRUE, 
                                  distmetric = Bruvo.distance, progress = TRUE)

#Getting alele diversity by population and alele
al_diversity <- alleleDiversity(data)

#Show results matrix
al_diversity$counts

# PCA analyses 

# PCA using simpleFreq (frecuencias alelicas) y Fst (distancias)

data <- deleteLoci(data, "F58red")
simple_freq <- simpleFreq(data, samples = Samples(data))
distances_Fst <- calcPopDiff(simple_freq, metric = "Fst")

pca <- cmdscale(distances_Fst, eig=T)

mycol <- c()
mycol[1:16] = "#f15544"
mycol[17:33] = "#6a93f3"
mycol[21] = "#77dd77"
mycol[7] = "black"
mycol[20] = "black"

apch <- c()
apch[1:16] = 16
apch[17:33] = 17

plot(-pca$points[,1], pca$points[,2],
     xlab = "PC1 (67.5%)",
     ylab = "PC2 (38.0%)",
     type = "p",
     bg = "black",
     pch = apch, #tipo marcador
     cex = 1.4, #tama?o puntos
     col = mycol
     #asignar color a grupos
)

legend("topleft", 
       legend = c("2x", "6x", "Kalksburg (2x)", "Kalksburg (6x)     ", "Kostalov"),
       pch = c(16, 17, 16, 17, 17),
       lty = NA,
       lwd = 2,
       cex = 0.8, 
       x.intersp = 0.4,
       col = c("#f15544", "#6a93f3", "black", "black", "#77dd77")
)

text(pca[,1], pca[,2], PopNames(data), cex = 0.6, pos = 4)


# PCA with genodive data (both cytotypes) Rho correlation matrix 

pca_gv <- read.table(file = "pca_gv.txt", header = TRUE)
head(pca_gv)

mycol <- c()
mycol[1:16] = "#f15544"
mycol[17:33] = "#6a93f3"
mycol[21] = "#77dd77"
mycol[7] = "black"
mycol[20] = "black"

apch <- c()
apch[1:16] = 16
apch[17:33] = 17

plot(pca_gv[,2], pca_gv[,3],
     xlab = "PC1 (53.6%)",
     ylab = "PC2 (6.9%)", 
     type = "p",
     bg = "black",
     pch = apch, 
     cex = 1.4,
     col = mycol
)

legend("topleft", 
       legend = c("2x", "6x", "Kalksburg (2x)", "Kalksburg (6x)     ", "Kostalov"),
       pch = c(16, 17, 16, 17, 17),
       lty = NA,
       lwd = 2,
       cex = 0.8, 
       x.intersp = 0.4,
       col = c("#f15544", "#6a93f3", "black", "black", "#77dd77")
)

#PCA mixed ploidy population with genodive data 

pca_gv <- read.table(file = "pca_gv_mixta.txt", header = TRUE)
head(pca_gv)

mycol <- c()
mycol[1:12] = "#f15544"
mycol[13:26] = "#6a93f3"

apch <- c()
apch[1:12] = 16
apch[17:26] = 17

plot(pca_gv[,2], pca_gv[,3],
     xlab = "PC1 (32.5%)",
     ylab = "PC2 (8,7%)", 
     type = "p",
     bg = "black",
     pch = apch, 
     cex = 1.4, 
     col = mycol
)

legend("topleft", 
       legend = c("Kalksburg (2x)", "Kalksburg (6x)      "),
       pch = c(16, 17),
       lty = NA,
       lwd = 2,
       cex = 0.8, 
       x.intersp = 0.4,
       col = c("#f15544", "#6a93f3")
)


