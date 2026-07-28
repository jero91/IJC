library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(FlowSorted.Blood.EPIC)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(ggplot2)
library(scales)
library(viridis)
library(readr)
library(tidyr)
library(gridExtra)
library(TCGAbiolinks)
library(biomaRt)
library(broom.mixed)
library(ggpubr)
library(plyranges)
#library(sqldf)
library(heatmap3)
library(dplyr)
library(gplots)
library(BioGeoBEARS)
library(tidyverse)

if(!require(installr)) {
  install.packages("installr"); require(installr)} 
updateR()  
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("FlowSorted.Blood.EPIC")
BiocManager::install("minfiData")
BiocManager::install("Gviz")
BiocManager::install("DMRcate")
BiocManager::install("readr")
BiocManager::install("rapport")
BiocManager::install("plyranges")
BiocManager::install("sqldf")
install.packages("dplyr")
install.packages("gplots")
install.packages('BioGeoBEARS')
install.packages('tidyverse')
install.packages("BiocManager")
install.packages("")


install.packages("readr")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")



CosmicSSEPIC <- CosmicCompleteCL[CosmicCompleteCL$GENE_NAME %in% "TCF4",]


getwd()

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
rownames(DMPs0.15Signi)

write.table(mSetSqFltAnnotated, file="mSetSqFltAnnotated.csv", sep=",", row.names=FALSE)

dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")

#list.files(dataDirectory, recursive = TRUE)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann850k[ann850k$Name %in% "cg09659072",]
ann450k[ann450k$Name %in% "cg09659072",]

ann850k <- as.data.frame(ann850k)
ann450k <- as.data.frame(ann450k)

colnames(ann850k)

ann850kSS <- ann850k[,c("Name","chr","Islands_Name","Relation_to_Island","UCSC_RefGene_Name", "UCSC_RefGene_Group", "GencodeBasicV12_NAME","GencodeBasicV12_Group" )]

#head(ann850k)

rm(targets)
#targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop.csv")


rgSet <- read.metharray.exp(base="/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/IdatsAutopsiasFinish34",  targets=targets)
#rgSet

# Save a single object to a file
saveRDS(rgSet, "rgSet.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DEL OBJETO RGSET EN RDS !!!!!!!!!!!!!
rgSet <- readRDS("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/rgSet.rds")

head(assay(rgSet))


detP <- detectionP(rgSet)
#head(detP)

pal <- brewer.pal(8,"Dark2")

#Para comparar la calidad de las muestras
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")


#Para filtrar por p-value para filtar muestras

keep <- colMeans(detP) < 0.01
rgSet <- rgSet[,keep]
rgSet
targets <- targets[keep,]
targets[,1:5]
detP <- detP[,keep]
dim(detP)
rgSet


# normalizacion 
mSetSq <- preprocessQuantile(rgSet) 

mSetRaw <- preprocessRaw(rgSet)

#Visualizar normalizacion
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))


#Otro metodo de normalizar:ssnoob
MSet.noob <- preprocessNoob(rgSet)
GRset.funnorm <- preprocessFunnorm(rgSet)


par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(MSet.noob), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))


#PCA por Sample_Group
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)



#Cambio sample-group por tissue y gender
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$TISSUE)])
legend("top", legend=levels(factor(targets$TISSUE)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$GENDER)])
legend("top", legend=levels(factor(targets$GENDER)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Rango.edades)])
legend("top", legend=levels(factor(targets$Rango.edades)), text.col=pal,
       bg="white", cex=0.7)
#Los porcentajes de los componentes principales son muy bajos.


#Cambiando los componentes que tenemos en cuenta:
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$GENDER)], dim=c(1,3))
legend("top", legend=levels(factor(targets$GENDER)), text.col=pal, 
       cex=0.7, bg="white")


# ensure probes are in the same order in the mSetSq and detP objects    
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
rowMeans(detP < 0.01) == ncol(mSetSq)

#Irene:
bad_pos <- row.names(as.data.frame(detP))[rowMeans(detP) > 0.01]
rgSet <- rgSet[minfi::featureNames(rgSet)[!(minfi::featureNames(rgSet) %in% bad_pos)], ]

mSetSq <- mSetSq[minfi::featureNames(mSetSq)[!(minfi::featureNames(mSetSq) %in% bad_pos)], ]
#Fin Irene

keep <- rowSums(detP < 0.01) == ncol(mSetSq) #Si el numero de columanas(muestras) que pasan el filtro (p-value < 0.01).
#Es decir me esta quitando todas las sondas que tienen al menos un valor no significativo(menor de 0.01).
mSetSq <- mSetSq[keep,]

table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
mSetSq

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in%   c("chrX","chrY")])
keep <- !(featureNames(mSetSq) %in% ann850k$Name[ann850k$chr %in%   c("chrX","chrY")])

table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- mSetSq[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

saveRDS(mSetSqFlt, "mSetSqFlt.rds") 
saveRDS(mSetSq, "mSetSq.rds")                    # !!!!!!!!!!Punto control!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!Punto control!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mSetSqFlt <- readRDS("mSetSqFlt.rds")
mSetSq <- readRDS("mSetSq.rds")

#Anotar genes
mSetSqFltAnnotated <- getAnnotation(mSetSqFlt)
head(mSetSqFltAnnotated)

write.table(mSetSqFltAnnotated, file="mSetSqFltAnnotated.csv", sep=",", row.names=FALSE)



# exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)


#PCA tras el filtrado
par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$TISSUE)])
legend("top", legend=levels(factor(targets$TISSUE)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",  
        col=pal[factor(targets$GENDER)])
legend("top", legend=levels(factor(targets$GENDER)), text.col=pal,
       bg="white", cex=0.7)


# calculate M-values for statistical analysis

mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

saveRDS(mVals, "mVals.rds")

mVals <- readRDS("mVals.rds")


mVals <- getM(mSetSq)


# calculate beta for statistical analysis

bVals <- getBeta(mSetSqFlt)
bVals <- getBeta(mSetSq)


head(bVals[,1:5])

# Save a single object to a file
saveRDS(bVals, "bVals.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DEL OBJETO bVals EN RDS. !!!!!!!!!!!!!
bVals <- readRDS("bVals.rds")


bValsDF <- as.data.frame(bVals)

# representamos m-value and betas
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$TISSUE, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$TISSUE)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$TISSUE, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$TISSUE)), 
       text.col=brewer.pal(8,"Dark2"))

summary(targets)

########                    LIMMA
#Creamos la matriz de diseno                            

# this is the factor of interest
tissue <- factor(targets$TISSUE)
# this is the individual effect that we need to account for 
gender <- factor(targets$GENDER) 
# use the above to create a design matrix
design <- model.matrix(~0+tissue+gender, data=targets) 
colnames(design) <- c(levels(tissue),levels(gender)[-1])
# fit the linear model 
#dim(mVals)
fit <- lmFit(mVals, design)
fit2 <- lmFit(bVals, design)
bVals
# create a contrast matrix for specific comparisons                     !!!esto no lo entiendo mucho y quiza sea impotante para los resultados.
contMatrix <- makeContrasts(CONTROL-COVID,
                            M-F,
                            levels=design) 
contMatrix


#Creamos la matriz de dise?o   sin genero porque si eliminamos las cromosomas X e Y no hay que ajustar por ese valor                         error dimension

# this is the factor of interest
tissue <- factor(targets$TISSUE)
# this is the individual effect that we need to account for 
#gender <- factor(targets$GENDER) 
# use the above to create a design matrix
design <- model.matrix(~0+tissue, data=targets) 
colnames(design) <- levels(tissue)
# fit the linear model 
#dim(mVals)
fit <- lmFit(bVals, design)
##fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons                     !!!esto no lo entiendo mucho y quiza sea impotante para los resultados.
contMatrix <- makeContrasts(CONTROL-COVID,
                            levels=design) 
contMatrix




fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))



#Este match es muy imprtante para que anote en orden
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann850kSub)
DMPsSigni <- subset(DMPs, adj.P.Val < 0.05)

#Coger de DMPs las que tengan p-value menor de 0.05 y entonces ejecuto mi escript de la media de la beta diferencia.
#Hago heatmap de las 4 clg que tengo y bajo a 0.15 y hago el heatmap tambien.

head(DMPs)
dim(DMPs)
saveRDS(DMPs, "DMPs.rds")
DMPs <- readRDS( "DMPs.rds")


DMPs1000 <- topTable(fit2, num=1000, coef=1,sort.by="logFC", genelist=ann850kSub)
head(DMPs1000)
DMPs100 <- topTable(fit2, num=100, coef=1,sort.by="logFC", genelist=ann850kSub)

write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

write.table(DMPs100, file="DMPs100.csv", sep=",", row.names=FALSE)

# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$TISSUE, ylab = "Beta values")
})

suma_renglones <- function (x) {
  y=integer(nrow(x))
  z=integer(nrow(x))
  r=integer(nrow(x))
  names <- colnames(x)
  for(i in 1:nrow(x)){
    for (j in 1:ncol(x)) {
      if(is.numeric(x[i,j]) && (names[j]=="204776850059_R01C01"
                                || names[j]=="204776850059_R05C01"
                                || names[j]=="204776850059_R07C01"
                                || names[j]=="204776850059_R08C01"
                                || names[j]=="204860850098_R05C01"
                                || names[j]=="204860850098_R07C01"
                                || names[j]=="204860850119_R01C01"
                                || names[j]=="204860850119_R02C01"
                                || names[j]=="204860850155_R01C01"
                                || names[j]=="204860850155_R04C01"
                                || names[j]=="204860850155_R06C01"
                                || names[j]=="204860850155_R07C01"
                                || names[j]=="205035150127_R05C01"
                                || names[j]=="205061410062_R03C01"
                                || names[j]=="205061410062_R07C01"
                                || names[j]=="205061410147_R01C01"
                                || names[j]=="205061410175_R03C01"
                                || names[j]=="205061410175_R07C01" )) 
      {
        y[i] = y[i]+(x[i,j])
      }else if (is.numeric(x[i,j]) && (names[j]=="204776850059_R06C01"
                                       || names[j]=="204860850098_R02C01"
                                       || names[j]=="204860850098_R06C01"
                                       || names[j]=="204860850098_R08C01"
                                       || names[j]=="204860850155_R02C01"
                                       || names[j]=="204860850155_R03C01"
                                       || names[j]=="204860850155_R05C01"
                                       || names[j]=="205035150127_R03C01"
                                       || names[j]=="205035150127_R07C01"
                                       || names[j]=="205061400131_R08C01"
                                       || names[j]=="205061410011_R01C01"
                                       || names[j]=="205061410062_R05C01"
                                       || names[j]=="205061410147_R03C01"
                                       || names[j]=="205061410175_R01C01"
                                       || names[j]=="205061410175_R05C01"
                                       || names[j]=="205061410062_R01C01")) {
        z[i] = z[i]+(x[i,j])
      }
    }
    dim(y)
    dim(x)
    y[i] <-y[i]/18
    z[i] <-z[i]/16
    r[i] = y[i]-z[i]
  }
  
  x = as.data.frame(x)
  x[, "mediaControl"] <- y
  x[, "mediaCovid"] <- z
  x[, "Diferencia"] <- r
  
  return(x)
}


BetaDiferenciaMedia = suma_renglones(bVals)

write.table(BetaDiferenciaMedia, file="BetaDiferenciaMedia.csv", sep=",", row.names=TRUE)

getwd()
# Save a single object to a file
saveRDS(BetaDiferenciaMedia, "BetaDiferenciaMedia.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DE. !!!!!!!!!!!!!
BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")

#BetaDiferenciaMedia[ rownames(BetaDiferenciaMedia) %in% "cg06055229",]


#BetaDiferenciaMediaSS <- subset(BetaDiferenciaMedia, Diferencia > 0.2 | Diferencia < -0.2)
#dim(BetaDiferenciaMediaSS)
#head(BetaDiferenciaMediaSS)

#head(BetaDiferenciaMediaSSAnot)


#mSetSqFltAnnotatedSS <- subset(mSetSqFltAnnotated, rownames(mSetSqFltAnnotated) == rownames(BetaDiferenciaMediaSS))
#mSetSqFltAnnotatedSS.2 <- mSetSqFltAnnotated[rownames(mSetSqFltAnnotated) %in% rownames(BetaDiferenciaMediaSS),]

#head(mSetSqFltAnnotated)
#head(mSetSqFltAnnotatedSS.2)

t <- 1
#for(i in  rownames(BetaDiferenciaMediaSS)){

#    print(t)
#    print(i)

#    mSetSqFltAnnotatedSS[t] <- subset(mSetSqFltAnnotated, rownames(mSetSqFltAnnotated) == i)

#   t <- t+1

#   }


#mSetSqFltAnnotatedSS <- subset(mSetSqFltAnnotated, rownames(mSetSqFltAnnotated) == "cg06295987")
#head(mSetSqFltAnnotatedSS)

#mSetSqFltAnnotated["cg06295987"]


#DMP.GUI(DMP=mSetSqFltAnnotatedSS,beta=bVals,pheno=DMPs$TISSUE)


#heatmap3::heatmap3(mSetSqFltAnnotatedSS)


BetaDiferenciaMediaSS <- subset(BetaDiferenciaMedia, Diferencia > 0.15 | Diferencia < -0.15)








#                       DMR                     #

?cpg.annotate
myAnnotationM <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "CONTROL - COVID", arraytype = "EPIC")

#Probar con beta-values
myAnnotationRS <- cpg.annotate(object = mSetSqFlt, datatype = "array",  
                               analysis.type = "differential", design = design, 
                               contrasts = TRUE, cont.matrix = contMatrix, 
                               coef = "CONTROL - COVID")

#Probar con beta-values
myAnnotationB <- cpg.annotate(object = bVals, datatype = "array",  what = "Beta", 
                              analysis.type = "differential", design = design, 
                              contrasts = TRUE, cont.matrix = contMatrix, 
                              coef = "CONTROL - COVID")

"CONTROL - COVID" %in% colnames(contMatrix)
colnames(contMatrix)

str(myAnnotation)


?dmrcate
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges


#Probando dmrcate
DMRs <- dmrcate(myAnnotationRS, lambda=1000, C=2,betacutoff = 0.1, min.cpgs = 5)




# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$TISSUE))]
names(groups) <- levels(factor(targets$TISSUE))
cols <- groups[as.character(factor(targets$TISSUE))]

# draw the plot for the top DMR
par(mfrow=c(2,2))
DMR.plot(ranges = results.ranges, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")
DMR.plot(ranges = results.ranges, dmr = 3, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")
DMR.plot(ranges = results.ranges, dmr = 5, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")
DMRs

# Save a single object to a file
saveRDS(DMRs, "DMRs.rds")
# Restore it under a different name     
DMRs <- readRDS("DMRs.rds")


results.ranges

?DMRcatedata
browseVignettes('DMRcatedata')







#####                       CHAMP                           #####

#BiocManager::install("ChAMP")
#BiocManager::install("DMP.GUI")
library(ChAMP)


DMP.GUI(DMP=DMPs[[1]],beta=myNorm,pheno=myLoad$pd$TISSUE)

DMP.GUI(DMP=DMPs,beta=bVals,pheno=DMPs$TISSUE)

CpG.GUI()


testDir <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopCHAMP"

testDir <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/IdatsAutopsiasFiltered30"

myLoad <- champ.load(testDir,arraytype="EPIC", method="minfi")

myLoad <- champ.load(testDir,arraytype="EPIC")

class(EPICSimData)

data(EPICSimData)
CpG.GUI(arraytype="EPIC")
champ.QC() # Alternatively QC.GUI(arraytype="EPIC")
myNorm <- champ.norm(arraytype="EPIC")
champ.SVD()
# If Batch detected, run champ.runCombat() here.This data is not suitable.
myDMP <- champ.DMP(arraytype="EPIC")
DMP.GUI()


?champ.import
?champ.load

list.files(testDir)
ss <- read.csv(paste0(testDir,"/SampleSheetAutop.csv"))
length(unique(ss$Basename))
length(ss$Basename)

champ.import(directory = testDir, offset=100, arraytype="EPIC")
?champ.filter




myLoad <- champ.load(testDir,arraytype="EPIC")

CpG.GUI(CpG=rownames(myLoad$beta),arraytype="EPIC")

head(myLoad$beta)

#Aqui sacar el data frame y poner los nombres. Voy a intentar cuando hago myNorm que pierdo los 
#nombres hacer que coja los nombres de pd$sampleID, luego filtro y me quedo con las cpg mas significativas y represento.

beta <- myLoad$beta
colnames(beta) <- myLoad$pd$SAMPLE.ID


myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=5)

# Save a single object to a file
saveRDS(myNorm, "myNorm.rds")
# Restore it under a different name             !!!!!!!COMIENZO !!!!!!!!!!!!!
myNorm <- readRDS("myNorm.rds")


myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$TISSUE,arraytype="EPIC")

DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$TISSUE)

myNorm <- champ.norm(beta=beta,arraytype="EPIC",cores=5)
myDMP2 <- champ.DMP(beta = myNormB,pheno=myLoad$pd$TISSUE,arraytype="EPIC")
DMP.GUI(DMP=myDMP2[[1]],beta=myNormB,pheno=myLoad$pd$TISSUE)






#Vamos a hacer un subset del objeto de champ(myDMP) para poder obtener los genes asociados y los heatmap

myDMP2SS <- myDMP$CONTROL_to_COVID[rownames(myDMP$CONTROL_to_COVID) %in% rownames(BetaDiferenciaMediaSS),]

head(myDMP2SS)

mSetSqFltAnnotatedSS.2 <- mSetSqFltAnnotated[rownames(mSetSqFltAnnotated) %in% rownames(BetaDiferenciaMediaSS),]




myDMP$CONTROL_to_COVID

# Save a single object to a file
saveRDS(myDMP2, "myDMP2.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DEL OBJETO RGSET EN RDS. !!!!!!!!!!!!!
myDMP <- readRDS("myDMP2.rds")

myLoad$pd$SAMPLE.ID

BiocManager::install("heatmap3")

heatMapMeth(myDMP)
heatmap3(myDMP$CONTROL_to_COVID$deltaBeta)

heatmap3(DMPs)
heatmap3(bVals)

head(myDMP)

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$TISSUE,method="Bumphunter",arraytype="EPIC")

myDMR <- champ.DMR(beta=myNormB,pheno=myLoad$pd$TISSUE,method="Bumphunter",arraytype="EPIC")

# Save a single object to a file
saveRDS(myDMR, "myDMR.rds")
# Restore it under a different name             !!!!!!!COMIENZO !!!!!!!!!!!!!
myDMR <- readRDS("myDMR.rds")


# Save a single object to a file
saveRDS(myLoad, "myLoad.rds")
# Restore it under a different name             !!!!!!!COMIENZO !!!!!!!!!!!!!
myLoad <- readRDS("myLoad.rds")


#Para ejecutar DMR.GUI necesito myNorm
DMR.GUI(DMR=myDMR,pheno=myLoad$pd$TISSUE)

matrix <- as.matrix(myDMP)

myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")

head(myGSEA$DMR)
head(myGSEA$DMP)


#Creo que no se puede hacer cna porque me esta cogienod una columna X que no existe.
myCNA <- champ.CNA(intensity=myLoad$intensity,pheno=myLoad$pd$TISSUE,arraytype="EPIC",controlGroup="CONTROL")

myCNA <- champ.CNA(intensity=as.data.frame(myLoad$intensity),pheno=myLoad$pd$TISSUE,arraytype="EPIC",controlGroup="CONTROL")


# Save a single object to a file
saveRDS(myCNA, "myCNA.rds")
# Restore it under a different name             !!!!!!!COMIENZO !!!!!!!!!!!!!
myCNA <- readRDS("myCNA.rds")

head(myCNA)
champ.SVD()

library(DNAcopy)


head(myCNA$groupResult)

head(myCNA$groupResult$COVID$seg.mean)


library(plot_segments)
plot_segments(
  myCNA,
  chromosomes = paste0("chr", c(1:22, "X", "Y")),
  max_Y_height = 6,
  circular = FALSE,
  cn = "absolute",
  highlight = myCNA$groupResult,
)

plot(myCNA$groupResult)


CNA(genomdat, chrom, maploc, data.type=c("logratio","binary"),
    sampleid=NULL, presorted = FALSE)

CNA.object <- CNA(myCNA$groupResult$COVID$seg.mean, myCNA$groupResult$COVID$chrom, myCNA$groupResult$COVID$loc.start)

smoothed.CNA.object <- smooth.CNA(CNA.object)

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)

segment.smoothed.CNA.object <- segment(smoothed.CNA.object)

plot(segment.smoothed.CNA.object, plot.type="w")

plot(segment.smoothed.CNA.object, plot.type="s")

"CNV.genomeplot"(object, chr = "all", chrX = TRUE, chrY = TRUE, centromere = TRUE, detail = TRUE, main = NULL, ylim = c(-1.25, 1.25), set_par = TRUE, cols = c("red", "red", "lightgrey", "green", "green"))


myCombat <- champ.runCombat(beta=myNorm,variablename = "TISSUE", pd=myLoad$pd,batchname=c("Slide"))

sampleResult <- myLoad$intensity
sampleResult[[colnames(intsqnlogratio)[i]]]
sampleResult[[colnames(intsqnlogratio)[i]]]
colnames(myLoad$pd)[19]
colnames(myLoad$pd)[19]
#https://rdrr.io/bioc/ChAMP/src/R/champ.CNA.R


intsqn <- normalize.quantiles(as.matrix(intensity))
intsqnlogratio <- apply(intsqnlog[,which(!pheno %in% controlGroup)],2,function(x) x - rowMeans(as.data.frame(intsqnlog[,which(pheno %in% controlGroup)])))



names <- colnames(myLoad$intensity)
names <- colnames(as.data.frame(myLoad$intensity))
names <- colnames(as.matrix(myLoad$intensity))

names


#########           conumee     ######
BiocManager::install("conumee")

library("conumee")
library("minfiData")
library("minfi")

RGsetTCGA <- read.EPIC.exp(base ="/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/IdatsAutopsiasFiltered30")  # use default parameters for vignette examples
read.

MsetTCGA <- preprocessIllumina(rgSet)

tcga.data <- CNV.load(MsetTCGA)

tcga.controls <- grep("-11A-", names(tcga.data))
names(tcga.data)


minfi.data <- CNV.load(MsetEx)
minfi.controls <- pData(MsetEx)$status == "normal"




#vAMOS HA HACER EL HEATMAP DE LA TAB?A CON 4 SONDAS:

bVals

bValsSS.2 <- bVals[rownames(bVals) %in% rownames(BetaDiferenciaMediaSS),]
head(bValsSS)

library(heatmap3)
heatmap3(bValsSS)

colnames(bValsSS) <- myLoad$pd$SAMPLE.ID

heatmap3(bValsSS.2)



#Ejecucion obtener los DMP con la informacion de los genes y los heatmap filtrando por betaDiferencia.

#Primero filtramos por valor de beta diferencia
BetaDiferenciaMedia0.2 <- subset(BetaDiferenciaMedia, Diferencia > 0.2 | Diferencia < -0.2)

#Aqui cogemos las betas para hacer el heatmap
bValsSS0.2 <- bVals[rownames(bVals) %in% rownames(BetaDiferenciaMedia0.2),]
dim(bValsSS0.19)
colnames(bValsSS0.2) <- myLoad$pd$SAMPLE.ID
heatmap3(bValsSS0.2)

heatmap3(bValsSS0.2, ColSideCut=10)



#Vamos a ver los nombres de los genes asociados s las sondas.
mSetSqFltAnnotatedSS0.15 <- mSetSqFltAnnotated[rownames(mSetSqFltAnnotated) %in% rownames(BetaDiferenciaMedia0.15),]
mSetSqFltAnnotatedSS0.15$UCSC_RefGene_Name
write.table(mSetSqFltAnnotatedSS0.15$UCSC_RefGene_Name, file="NombreGenes.csv", sep=",", row.names=TRUE)

getwd()

#Nueva ejecucion con 34 muestras:
# Save a single object to a file
saveRDS(BetaDiferenciaMedia, "BetaDiferenciaMedia.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DE. !!!!!!!!!!!!!
BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")

DMPs <- readRDS( "DMPs.rds")

getwd()



#Filtrado sondas significativas:


BetaDiferenciaMedia0.15 <- subset(BetaDiferenciaMedia, Diferencia > 0.15 | Diferencia < -0.15)
dim(BetaDiferenciaMedia0.15)
rownames(DMPs)
dim(DMPs)

DMPs
head(DMPs)

DMPs0.15 <- DMPs[match(rownames(BetaDiferenciaMedia0.15),rownames(DMPs)),]

#Como funciona el match: El primer nombre de betasDiferencia esta en la posicion 16 de DMPs,
#entonces podemos eldataframe mas grande al principio y al final.
#BetaDiferenciaMedia0.15[1,]
#DMPs[16,]


DMPs0.15[DMPs0.15$Name %in% "cg27496650", ]
DMPs0.15[DMPs0.15$Name %in% "cg08884752", ]

DMPs0.15Signi <- subset(DMPs0.15, adj.P.Val < 0.05) #Los p.valores son difentees antes y despues del filtrado!!!! no usar subset!

dim(DMPs0.15Signi)

DMPs0.15Signi[DMPs0.15Signi$Name %in% "cg27496650", ]
DMPs0.15Signi[DMPs0.15Signi$Name %in% "cg08884752", ]



DMPs0.15Signi <- DMPs0.15[DMPs0.15$adj.P.Val < 0.05,]

BetaDiferenciaMedia0.15Signi <- BetaDiferenciaMedia0.15[rownames(BetaDiferenciaMedia0.15) %in% rownames(DMPs0.15Signi),]

BetaDiferenciaMedia0.15Signi <- BetaDiferenciaMedia0.15Signi[,c("mediaControl","mediaCovid")]
write.table(BetaDiferenciaMedia0.15Signi, file="BetaDiferenciaMedia0.15Signi.csv", sep=",", row.names=TRUE)


#match(rownames(DMPsSigni),rownames(BetaDiferenciaMedia0.15))


BetaDiferenciaMedia0.15

bVals0.15Signi <- bVals[match(rownames(DMPs0.15Signi),rownames(bVals)),]

bVals0.15Signi <- na.omit(bVals0.15Signi)

dim(bVals0.15Signi)
colnames(bVals0.15Signi) <- targets$EPIC.Sample.Name
heatmap3(bVals0.15Signi, ColSideCut=10)
heatmap3(bVals0.15Signi)

targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")



BetaDiferenciaMedia0.15Signi <- 
  
  #filtrado2:
  
  #DMPsSigni <- subset(DMPs, adj.P.Val < 0.05)
  DMPsSigni <- DMPs[DMPs$adj.P.Val < 0.05,]

DMPsSigni0.15 <- DMPsSigni[match(rownames(BetaDiferenciaMediaSS),rownames(DMPsSigni)),]

DMPsSigni0.15 <- DMPsSigni0.15[!is.na(DMPsSigni0.15$chr),]

match(rownames(DMPsSigni0.15),rownames(DMPs0.15Signi))  







#No hace falta??
DMPs0.15Signi <- DMPs[match(rownames(bVals0.15Signi),rownames(DMPs)),]
dim(DMPs0.15)
listaGenes <- DMPs0.15$UCSC_RefGene_Name
write.table(listaGenes, file="listaGenes.csv", sep=",", row.names=TRUE)





#Irene:

IreneProbes <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/IreneResultsProbes.csv", TRUE, sep = ",")
dim(IreneProbes)

IreneProbes$Name
rownames(DMPs0.15Signi)

JOIN <- IreneProbes[match(IreneProbes$Name,rownames(DMPs0.15Signi)),]

JOINnotNA <- na.omit(JOIN)
dim(JOINnotNA)


#Calcular cuantas DMPs hay significativo.
DMPsSignif <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPsSignif)

#Compara irene
DMPs <- readRDS( "DMPs.rds")
setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

# Save a single object to a file
saveRDS(DMPs0.15Signi, "DMPs0.15Signi.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DE. !!!!!!!!!!!!!
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
write.table(DMPs0.15Signi, file="DMPs0.15Signi.csv", sep=",", row.names=TRUE)

heatmap3(DMPs0.15Signi)


nombreGenesUSDC <- DMPs0.15Signi$UCSC_RefGene_Name
nombreGenesUSDC <- na.omit(nombreGenesUSDC)
head(nombreGenesUSDC)
write.table(DMPs0.15Signi$UCSC_RefGene_Name, file="listaGenes.csv", sep=",", row.names=TRUE)

nombreGenesUSDCSplit<- strsplit(as.character(nombreGenesUSDC), ";")

#nombreGenesUSDCSplitUnicos <- unique(nombreGenesUSDCSplit)

nombreGenesUSDCSplitUnicosUnlist <-unlist(nombreGenesUSDCSplit)

nombreGenesUSDCSplitUnicos <- unique(nombreGenesUSDCSplitUnicosUnlist)

nombreGenesUSDCSplitUnicos <- as.data.frame(nombreGenesUSDCSplitUnicos)
write.table(nombreGenesUSDCSplitUnicos, file="listaGenesUSDC.csv", sep=",", row.names=TRUE)

nombreGenesGencode <- DMPs0.15Signi$GencodeBasicV12_NAME

nombreGenesGencodeSplit<- strsplit(as.character(nombreGenesGencode), ";")

nombreGenesGencodeSplitUnlist <-unlist(nombreGenesGencodeSplit)

nombreGenesGencodeSplitUnlistUnicos <- unique(nombreGenesGencodeSplitUnlist)

nombreGenesGencode <- as.data.frame(nombreGenesGencodeSplitUnlistUnicos)

write.table(nombreGenesGencode, file="listaGenesGencode.csv", sep=",", row.names=TRUE)

#INNER JOIN
genesJoin <- nombreGenesUSDCSplitUnicos[nombreGenesUSDCSplitUnicos  %in% nombreGenesGencode,]


#Isla o promotor:
relationToIsland <- unique(DMPs0.15Signi$Relation_to_Island)

nombresIslas <- c("N_Shore","OpenSea")

cuantasIslas <- DMPs0.15Signi$Relation_to_Island[DMPs0.15Signi$Relation_to_Island %in% nombresIslas]

#cuantasIslas <- count([DMPs0.15Signi$Relation_to_Island %in% nombresIslas,])

#radar plot
plot <- ggplot(as.data.frame(DMPs0.15Signi$Relation_to_Island)) + geom_bar()

plot <- ggplot(DMPs0.15Signi, aes(Relation_to_Island, frequency(Relation_to_Island), fill = Relation_to_Island)) +
  geom_bar(width = 1, stat = "identity", color = "white") 



plot <- ggplot(DMPs0.15Signi, aes(Relation_to_Island, frequency(Relation_to_Island), fill = Relation_to_Island)) +
  geom_bar(aes(y=count(Relation_to_Island)/120)) 


plot <- ggplot(DMPs0.15Signi, aes(Relation_to_Island, frequency(Relation_to_Island), fill = Relation_to_Island)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  geom_text(aes( label = scales::percent(frequency(Relation_to_Island)),
                 y=  frequency(Relation_to_Island)), stat= "identity", vjust = -.5)


plot
plot + coord_polar()

#Lluis, quiza con table


plot <-ggplot(DMPs0.15Signi, aes(x= Relation_to_Island)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..), y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="AGE") +
  scale_y_continuous(labels=percent)


plot <-ggplot(DMPs0.15Signi, aes(x=as.numeric(Relation_to_Island), y = ..prop..)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..), y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="AGE") +
  scale_y_continuous(labels=percent)
plot

grafico <- ggplot(DMPs0.15Signi, aes(x=as.numeric(Relation_to_Island), colour = rainbow(1))) + 
  geom_bar(aes(y = ..prop.. , fill = as.numeric(Relation_to_Island)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="Relation_to_Island") +
  scale_y_continuous(labels=percent)


grafico <- ggplot(DMPs0.15Signi, aes(x=Relation_to_Island)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = Relation_to_Island), stat="count", color = "white") +
  geom_text(aes( label = scales::percent((..count..)/sum(..count..)),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="Relation_to_Island") +
  scale_y_continuous(labels=percent)

grafico + coord_polar()
grafico
warning
plot <-ggplot(DMPs0.15Signi, aes(x= Relation_to_Island,  y = P.Value)) + geom_col()
plot


plot <-ggplot(DMPs0.15Signi, aes(x= Relation_to_Island) + 
                geom_bar(aes(y= (..count..)/sum(..count..))) +   scale_y_continuous(labels=percent_format()))

plot
plot <- ggplot(DMPs0.15Signi, aes(x= as.numeric(Relation_to_Island), y=..count..))
plot                

plot <- ggplot(data=DMPs0.15Signi, aes(x=as.numeric(Relation_to_Island))) + geom_line(aes(fill=..count..), stat="count", binwidth=10)



grafico <- ggplot(DMPs0.15Signi, aes(x=UCSC_RefGene_Group)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = Relation_to_Island), stat="count", color = "white") +
  geom_text(aes( label = scales::percent((..count..)/sum(..count..)),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="UCSC_RefGene_Group") +
  scale_y_continuous(labels=percent)
grafico + coord_polar()


grafico <- ggplot(DMPs0.15Signi, aes(x=Relation_to_Island)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = Relation_to_Island), stat="count", color = "white") +
  geom_text(aes(label = scales::percent((..count..)/sum(..count..)),
                y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="Relation_to_Island") +
  scale_y_continuous(labels=percent)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

grafico + coord_polar()

grafico <- ggplot(DMPs0.15Signi, aes(x=GencodeBasicV12_Group)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = Relation_to_Island), stat="count", color = "white") +
  geom_text(aes( label = scales::percent((..count..)/sum(..count..)),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="GencodeBasicV12_Group") +
  scale_y_continuous(labels=percent)
grafico + coord_polar()


grafico <- ggplot(DMPs0.15Signi, aes(x=Relation_to_Island)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = Relation_to_Island), stat="count", color = "white") +
  geom_text(aes(label = scales::percent((..count..)/sum(..count..), colour = "grey"),
                y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="Relation_to_Island") +
  scale_y_continuous(labels=percent)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "black", size = 0.5), 
        panel.grid = element_line(colour = "black", size = 5), 
        panel.grid.minor = element_line(colour = "black", size = 20))

grafico + coord_polar()




GrupoGenesUSDCSplit<- strsplit(as.character(DMPs0.15Signi$UCSC_RefGene_Group), ";")
GrupoGenesUSDCSplitUnlist <-unlist(GrupoGenesUSDCSplit)
GrupoGenesUSDCSplitUnlistDF <- as.data.frame(GrupoGenesUSDCSplitUnlist)

GrupoGenesUSDCUnique <- unique(DMPs0.15Signi$UCSC_RefGene_Group)

DMPs0.15Signi$UCSC_RefGene_Group[3]

result <- 0
rm(val)
count=1
for(val in DMPs0.15Signi$UCSC_RefGene_Group){
  
  GrupoGenesUSDCSplit<- strsplit(as.character(val), ";")
  GrupoGenesUSDCSplitUnlist <- unlist(GrupoGenesUSDCSplit[1])
  result[count] <- GrupoGenesUSDCSplitUnlist[1]
  count= count + 1
}

result <- na.omit(result)
resultDF <- as.data.frame(result)

GrupoGenesUSDCSplit<- strsplit(as.character(DMPs0.15Signi$UCSC_RefGene_Group[3]), ";")
val[2] <- GrupoGenesUSDCSplit[1]
GrupoGenesUSDCSplitUnlist <- unlist(GrupoGenesUSDCSplit[1])
GrupoGenesUSDCSplitUnlist[1]



grafico <- ggplot(GrupoGenesUSDCSplitUnlistDF, aes(x=GrupoGenesUSDCSplitUnlist)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = GrupoGenesUSDCSplitUnlist), stat="count", color = "white") +
  geom_text(aes( label = scales::percent((..count..)/sum(..count..)),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="GrupoGenesUSDCSplitUnlist") +
  scale_y_continuous(labels=percent)
grafico + coord_polar()

listaGenesGencode <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/listaGenesGencode.csv")

listaGenesUSDC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/listaGenesUSDC.csv")

genesJoin <- listaGenesUSDC[listaGenesUSDC$nombreGenesUSDCSplitUnicos %in% listaGenesGencode$nombreGenesGencodeSplitUnlistUnicos, ]

write.table(genesJoin, file="genesJoin.csv", sep=",", row.names=TRUE)

grafico <- ggplot(DMPs0.15Signi, aes(x=Relation_to_Island)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = Relation_to_Island), stat="count", color = "white", width = 0.75) +
  geom_text(aes(label = scales::percent((..count..)/sum(..count..), colour = "grey"),
                y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="Relation_to_Island") +
  scale_y_continuous(labels=percent)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "black", size = 0.5), 
        panel.grid = element_line(colour = "black", size = 5), 
        panel.grid.minor = element_line(colour = "black", size = 20))

grafico + coord_polar()


grafico <- ggplot(resultDF, aes(x=result)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = result), stat="count", color = "white", width = 0.75) +
  geom_text(aes(label = scales::percent((..count..)/sum(..count..), colour = "grey"),
                y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="result") +
  scale_y_continuous(labels=percent)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "black", size = 0.5), 
        panel.grid = element_line(colour = "black", size = 5), 
        panel.grid.minor = element_line(colour = "black", size = 20))

grafico + coord_polar()


grafico <- ggplot(DMPs0.15Signi, aes(x=Relation_to_Island)) + 
  geom_bar(aes(y = (..count..)/sum(..count..) , fill = Relation_to_Island), stat="count", color = "black", width = 0.75, lwd = 0.71) +
  geom_text(aes(label = scales::percent((..count..)/sum(..count..)),
                y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="Relation_to_Island") +
  scale_y_continuous(labels=percent)+
  scale_fill_viridis(discrete = TRUE, option = "A") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "black", size = 0.5))


grafico + coord_polar()


#Resultados GSEA

ResultsGSEA <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/resultGSEA.csv", TRUE, sep = ",")
ResultsGSEA <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/ResultsGSEAold.csv", TRUE, sep = ",")
colnames(ResultsGSEA)
plot <-ggplot(ResultsGSEA, aes(x=Gene.Set.Name,  y = X..Genes.in.Gene.Set..K.))
plot


barplot_GSEA <- function(data = data, title = "") {
  data$Gene.Set.Name <- gsub("^.*?_","", data$Gene.Set.Name)
  data$Gene.Set.Name <- factor(data$Gene.Set.Name, levels = data$Gene.Set.Name[order(as.numeric(data$X..Genes.in.Gene.Set..K.), decreasing = FALSE)])
  ggplot(data, aes_string(x = "X..Genes.in.Overlap..k.", y = "Gene.Set.Name", fill = "p.value")) +
    geom_col()+
    theme_bw() +
    scale_fill_continuous(low="red", high="red4", name = "p.value",
                          guide=guide_colorbar(reverse=TRUE)) +
    theme(axis.text = element_text(size = 12)) +
    ggtitle(title) +
    xlab(NULL) + ylab(NULL)
}

library(ggplot2)
plot <- barplot_GSEA(ResultsGSEA)
plot


#TCGA biolinks 
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
BiocManager::install("datatable")
library(datatable)
library(DT)

#Lanzamos la query a TCGA:
GDCqueryMetilacion <- GDCquery(project="TCGA-LUSC", data.category="DNA Methylation", sample.type= "Solid Tissue Normal")

GDCqueryExpresion <- GDCquery(project="TCGA-LUSC", data.category="Transcriptome Profiling", data.type="Gene Expression Quantification", sample.type= "Solid Tissue Normal")


#write.table(DMPs0.15Signi$Name, file="sondasSignificativas.csv", sep=",", row.names=FALSE)

##datatable(getResults(GDCqueryMetilacion), 
#          filter = 'top',
#          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#          rownames = FALSE)

common.patients <- intersect(substr(getResults(GDCqueryMetilacion, cols = "cases"), 1, 12),
                             substr(getResults(GDCqueryExpresion, cols = "cases"), 1, 12))

#Comprobar que funciona tambien con cases.submitter id.

query.met <- GDCquery(project = "TCGA-LUSC", data.category = "DNA Methylation",
                      sample.type= "Solid Tissue Normal" , barcode = common.patients[1:7])

query.exp <- GDCquery(project = "TCGA-LUSC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      sample.type= "Solid Tissue Normal",
                      barcode = common.patients[1:7])




datatable(getResults(query.met, cols = c("data_type","cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

datatable(getResults(query.exp, cols = c("data_type","cases")), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

data.met <- GDCdownload(query.met)


#dataMetilacion <- read.table(file = "jhu-usc.edu_LUSC.HumanMethylation450.3.lvl-3.TCGA-22-5472-11A-11D-1633-05.gdc_hg38.txt", header = TRUE, sep='\t')
#Hariamos 7 lineas como la anterior, una por cada muestra.


data.exp <- GDCdownload(query.exp)

#query.exp <- query.exp[[1]][[1]]

query.exp <- subset(query.exp[[1]][[1]], analysis_workflow_type == "HTSeq - FPKM-UQ")


dataExpresion <- read.table(file = "cd961e28-af13-4780-a068-35ff5548a3d2.FPKM-UQ.txt", header = TRUE, sep='\t')
#biolink leer muchos ficheros

#Creamos la tabla con los datos de metilacion, expresion y el nombre de los genes por sonda

BiocManager::install("biomaRt")
library(biomaRt)
#hgnc_swissprot <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','uniprotswissprot'),filters = 'ensembl_gene_id', values = 'ENSG00000139618', mart = ensembl)
#hgnc_swissprot

ensembl = useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
#ensembl = useMart(biomart = "ensembl", dataset= "mmusculus_gene_ensembl")


charg2 <- sapply(strsplit(dataExpresion$ENSG00000242268.2, '.', fixed=T), function(x) x[1])

dataExpresion$ENSG00000242268.2 <- sapply(strsplit(dataExpresion$ENSG00000242268.2, '.', fixed=T), function(x) x[1])

nombreGenes <- getBM(attributes='hgnc_symbol', 
                     filters = 'ensembl_gene_id', 
                     values = charg2, 
                     mart = ensembl)

nombreGenes <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), 
                     filters = 'ensembl_gene_id', 
                     values = charg2, 
                     mart = ensembl)



library(org.Hs.eg.db)

annots <- select(org.Hs.eg.db, keys=charg2, 
                 columns="SYMBOL", keytype="ENSEMBL")



#Hay 3 id que estan repetidos, elimino los 3 restantes, quiza seria mejor desdoblarlos o ver cuales son y porque!!!!
resdata2 <- resdata2[!duplicated(resdata2$ensembl_gene_id),]






#Como dataMetilacion tiene muchos genes por sonda me quedo con el primero, aunque seria mejor desdoblar.
#MetilacionUnlist <- unlist(data1$Gene_Symbol)

Gene_Symbol<- strsplit(as.character(dataMetilacion$Gene_Symbol), ";")
#Gene_SymbolDF <- as.data.frame(Gene_Symbol)

listaGenes <- c()
count <- 1
for(i in Gene_Symbol){
  listaGenes[count] <- i[1]
  count= count + 1
}

dataMetilacion$Gene_Symbol <- listaGenes



#dataMetilacion$expresion <- dataExpresion[dataExpresion$genes %in% dataMetilacion$Gene_Symbol, "X1095.57324501"]

#detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

#match(data1$Gene_Symbol,resdata2$hgnc_symbol)

#Con la tabla nombregenes y dataExpresion mapeamos por id del trascrito.
dataExpresion$hgnc_symbol <- nombreGenes[match(dataExpresion$V1,nombreGenes$ensembl_gene_id), "hgnc_symbol"]


dataMetilacion$expresion <- dataExpresion[match(dataMetilacion$Gene_Symbol,dataExpresion$hgnc_symbol), "V2"]

saveRDS(data1, "data1.rds")






#Cuantas de mis 120 sondas hay en los experimentos descargamos de TCGA(450K).


DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
#match(rownames(DMPs0.15Signi),data1$Composite.Element.REF)
#rownames(DMPs0.15Signi) %in% data1$Composite.Element.REF


result <- dataMetilacion[match(rownames(DMPs0.15Signi),dataMetilacion$Composite.Element.REF),]

result$promotor <- DMPs0.15Signi[match(rownames(DMPs0.15Signi),result$Composite.Element.REF),"UCSC_RefGene_Group"]
#DMPs0.15Signi$UCSC_RefGene_Group


#result$promotor <- DMPs0.15Signi[match(dataMetilacion$Composite.Element.REF,rownames(DMPs0.15Signi)),"UCSC_RefGene_Group"]

#result$promotor <- DMPs0.15Signi[rownames(DMPs0.15Signi) %in% data1$Composite.Element.REF, "UCSC_RefGene_Group"]


#Representamos la metliacion vs la expresion para las 50 sondas que tenemos informacion.


ggplot(result,aes(result$Beta_value, result$expresion)) +
  geom_point() +
  geom_smooth(method='lm')



#Hacemos el logaritmo de la expresion

result$expresion <- log(result$expresion)

#Separamos las que estan en promotor y no

resultNoPromotores <- subset(result,  promotor == "Body" | promotor == "3'UTR" )

resultPromotores <- subset(result,  promotor == "5'UTR" | promotor == "TSS1500" | promotor == "TSS200" | promotor == "1stExon")




ggplot(resultPromotores,aes(resultPromotores$Beta_value, resultPromotores$expresion)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(resultNoPromotores,aes(resultNoPromotores$Beta_value, resultNoPromotores$expresion)) +
  geom_point() +
  geom_smooth(method='lm')


#Representamos las 7 muestras para cada cpg.
#Empezamos por la sonda de CALCB. cg16394018

cg16394018.2 <- subset(dataMetilacion, Composite.Element.REF=="cg16394018" , select = c(Composite.Element.REF,Beta_value,expresion))

cg16394018 <- rbind(cg16394018,cg16394018.2)

common.patients[2]
#TCGA-22-5478
dataMetilacion <- read.table(file = "datosMetilacionTCGA/jhu-usc.edu_LUSC.HumanMethylation450.3.lvl-3.TCGA-22-5478-11A-11D-1633-05.gdc_hg38.txt", header = TRUE, sep='\t')

dataExpresion <- read.table(file = "datosExpresionTCGA/8ff70f54-8735-40b1-a556-b96e8889cc0f.FPKM-UQ.txt", header = FALSE, sep='\t')


dataExpresion$V1 <- sapply(strsplit(dataExpresion$V1, '.', fixed=T), function(x) x[1])


ensembl = useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

nombreGenes <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), 
                     filters = 'ensembl_gene_id', 
                     values = dataExpresion$V1, 
                     mart = ensembl)

write.table(cg16394018, file="cg16394018.csv", sep=",", row.names=FALSE)


#Automatizar:
common.patients[3]
#TCGA-22-5482
nombreArchivoMetilacion <- dataMetilacion[match(common.patients[3],dataMetilacion$Composite.Element.REF),"file_name"]

nombreArchivoMetilacion <- dataMetilacion[common.patients[3] %in% GDCqueryMetilacion$cases, "file_name"]

#GDCqueryMetilacion$results$cases ¿?

data <- GDCprepare(GDCqueryMetilacion)
data.met <- GDCdownload(query.met)
data.met.prepare <- GDCprepare(query.met)

query.met <- query.met[[1]][[1]]

query.exp <- query.exp[[1]][[1]]

cg16394018 <- read.csv("cg16394018.csv")

nombreArchivoMetilacion <- query.met[match(common.patients[3],query.met$cases.submitter_id),"file_name"]


#1 Query metilacion y expresion.
#2 Pacientes(muestras) comunes
#3 Query y descarga de datos de pacientes comunes
#4 Creamos la tabla con los nombres de los genes.
#5 Para cada uno de los pacientes comunes:
#6 obtenemos los nombres del archivos
#7 cargamos y tratamiento a los datos de expresion, le quitamos el sobrante del id.
#8 nos quedamos solo con el primer gen de la tabla de metilacion y primer promotor
#9 Mapeamos por nombres de genes y añadimos la expresion y el identificador de  muestra.
#10 Unimos todo en una tabla por cpg

cg22980156 <- c()

#i <- common.patients[1]

for(i in common.patients){
  
  nombreArchivoMetilacion <- query.met[match(i,query.met$cases.submitter_id),"file_name"]
  nombreArchivoExpresion <- query.exp[match(i,query.exp$cases.submitter_id),"file_name"]
  
  rutaMetilacion <- paste("datosMetilacionTCGA/",nombreArchivoMetilacion , sep = "")
  
  rutaExpresion <- paste("datosExpresionTCGA/",nombreArchivoExpresion , sep = "")
  
  dataMetilacion <- read.table(file = rutaMetilacion, header = TRUE, sep='\t')
  
  dataExpresion <- read.table(file = rutaExpresion, header = FALSE, sep='\t')
  
  dataExpresion$V1 <- sapply(strsplit(dataExpresion$V1, '.', fixed=T), function(x) x[1])
  
  
  Gene_Symbol<- strsplit(as.character(dataMetilacion$Gene_Symbol), ";")
  #Gene_SymbolDF <- as.data.frame(Gene_Symbol)
  
  listaGenes <- c()
  count <- 1
  for(j in Gene_Symbol){
    listaGenes[count] <- j[1]
    count= count + 1
  }
  
  dataMetilacion$Gene_Symbol <- listaGenes
  
  dataExpresion$hgnc_symbol <- nombreGenes[match(dataExpresion$V1,nombreGenes$ensembl_gene_id), "hgnc_symbol"]
  
  dataMetilacion$expresion <- dataExpresion[match(dataMetilacion$Gene_Symbol,dataExpresion$hgnc_symbol), "V2"]
  
  result <- dataMetilacion[match(rownames(DMPs0.15Signi),dataMetilacion$Composite.Element.REF),]
  
  result$promotor <- DMPs0.15Signi[match(rownames(DMPs0.15Signi),result$Composite.Element.REF),"UCSC_RefGene_Group"]
  
  result$promotor <- sapply(strsplit(as.character(result$promotor), ";"), function(x) x[1])
  
  cg22980156.2 <- subset(dataMetilacion, Composite.Element.REF=="cg22980156" , select = c(Composite.Element.REF,Beta_value,expresion))
  
  cg22980156.2$cases.submitter_id <- i
  
  cg22980156 <- rbind(cg22980156,cg22980156.2)
  
}

cg22980156$expresion <- log(cg22980156$expresion)

ggplot(cg22980156,aes(cg22980156$Beta_value, cg22980156$expresion)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(cg22980156,aes(cg22980156$Beta_value, cg22980156$expresion)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_fit_glance(method = 'lm')




#Inicio script:
#Lanzamos la query a TCGA:
GDCqueryMetilacion <- GDCquery(project="TCGA-LUSC", data.category="DNA Methylation", sample.type= "Solid Tissue Normal")

GDCqueryExpresion <- GDCquery(project="TCGA-LUSC", data.category="Transcriptome Profiling", data.type="Gene Expression Quantification", sample.type= "Solid Tissue Normal")

common.patients <- intersect(substr(getResults(GDCqueryMetilacion, cols = "cases"), 1, 12),
                             substr(getResults(GDCqueryExpresion, cols = "cases"), 1, 12))

#Comprobar que funciona tambien con cases.submitter id.

query.met <- GDCquery(project = "TCGA-LUSC", data.category = "DNA Methylation",
                      sample.type= "Solid Tissue Normal" , barcode = common.patients[1:7])

query.exp <- GDCquery(project = "TCGA-LUSC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      sample.type= "Solid Tissue Normal",
                      barcode = common.patients[1:7])


data.met <- GDCdownload(query.met)

data.exp <- GDCdownload(query.exp)

query.exp <- subset(query.exp[[1]][[1]], analysis_workflow_type == "HTSeq - FPKM-UQ")

query.met <- query.met[[1]][[1]]





#Sacamos el nombre de las 120 sondas significativas para introducirlas en methylationDB.

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
sondasCS <- c()
for(i in DMPs0.15Signi$Name){
  sondasCS <- paste(sondasCS,i, sep = "")
  sondasCS <- paste(sondasCS,",", sep = "")
}
sondasCS

rownames(DMPs0.15Signi)

betas_result <- read.csv("betas_result_2021-08-06-10-10.csv")

betas_resultSS <- read.csv("betas_result_2021-08-06-10-10.csv", skip=7)  #Epic de pulmon


betas_result450K <- read.csv("betas_result_2021-08-06-13-40.csv")

betas_result450KSS <- read.csv("betas_result_2021-08-06-13-40.csv", skip=7) #450K de pulmon


#Sript Lluis abrir archivos y ponerlos en otro directorio:

for(x in query.exp$id){
  new_directory_exp <- paste("/Users/lgimenez/Desktop/GDCdata/TCGA-LUSC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/", x, sep = "")
  setwd(new_directory_exp)
  list.of.files <- list.files(new_directory_exp,)
  file.copy(list.of.files, "/Users/lgimenez/Documents/datosExpresionTCGA")
}





#Cosmic:


listaGenesUSDC <- read.csv("listaGenesUSDC.csv")

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

unique(CosmicCompleteCL$SAMPLE_NAME)


betas_result450KSSUnique <- betas_result450KSS[!duplicated(betas_result450KSS$Sample.Name),]

CosmicSS450K <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_result450KSS$Sample.Name,]

#Investigar DIFERENCIAA!!! con %in%  ¿Duplicados?
#CosmicSS450KMerged <- merge(CosmicCompleteCL,betas_result450KSS, by.x = "SAMPLE_NAME", by.y = "Sample.Name" )

#duplicated(CosmicSS450KMerged$SAMPLE_NAME)

#CosmicSS450KMergedUnique <-CosmicSS450KMerged[!duplicated(CosmicSS450KMerged),]


CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% listaGenesUSDC$X,]


CosmicCTNND1 <- CosmicCompleteCL[CosmicCompleteCL$GENE_NAME %in% "CTNND1",]



#Para añadir los datos de metilacion a los datos de cosmic
#Extraemos los nombres d elos genes, y hacemos el unique.
#Toncees para cada uno de los unique que haga un for que rellene un subset con los datos que coincidan para ese gen.

CosmicSS450KNombreGenes <- CosmicSS450KGenes$GENE_NAME
CosmicSS450KNombreGenes <- unique(CosmicSS450KNombreGenes)

CosmicSS450KNombreGenes <- CosmicSS450K$GENE_NAME
CosmicSS450KNombreGenes <- unique(CosmicSS450KNombreGenes)
#j <- c()
#tablasGenes <- c()
#tablasGenesDF <- as.data.frame(tablasGenes)
#count <- 1
for(i in CosmicSS450KNombreGenes){
  
  #tablasGenes2 <- CosmicSS450KGenes[CosmicSS450KGenes$GENE_NAME %in% i,]
  
  #tablasGenes <- list(tablasGenes,tablasGenes2)
  
  assign(paste0("gen", i), CosmicSS450KGenes[CosmicSS450KGenes$GENE_NAME %in% i,])
  
  # tablasGenes[count] <- CosmicSS450KGenes[match(CosmicSS450KGenes$GENE_NAME,i),]
  
  #count <- count + 1
}
#tablasGenes <- na.omit(tablasGenes)
#tablasGenes[[1]]
#tablasGenesDF <- as.data.frame(tablasGenes)


#Tenemos datos de expresion y metialcion por genes, vamos a unir los datos por muestras:
#genSKI$expresion <- CosmicSS450KGenes[CosmicSS450KGenes$SAMPLE_NAME  %in% betas_result450KSS$Sample.Name ,"cg08884752"]   


b <- betas_result450KSS[betas_result450KSS$Sample.Name %in% genSKI$SAMPLE_NAME ,"cg08884752"]
a <- betas_result450KSS[match(CosmicSS450KGenes$SAMPLE_NAME,nombres),"cg08884752"]
#genSKI$expresion <- betas_result450KSS[match(CosmicSS450KGenes$SAMPLE_NAME,betas_result450KSS$Sample.Name),"cg08884752"]


rm(genSKI)

betas_result450KSSUnique <- betas_result450KSS[!duplicated(betas_result450KSS$Sample.Name),]

#a <- betas_result450KSS[match(betas_result450KSS$Sample.Name,genSKI$SAMPLE_NAME ),"cg08884752"]

expresionGen <- betas_result450KSSUnique[betas_result450KSSUnique$Sample.Name %in% genSKI$SAMPLE_NAME ,c("Sample.Name", "cg08884752")]
finalSKI <- merge(genSKI, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")



a <- betas_result450KSS[match(nombres,genSKI$SAMPLE_NAME ),"cg08884752"]
a <- na.omit(a)
genSKI$expresion <- a



library(TCGAbiolinks)
library(ggplot2)
library(biomaRt)
library(broom.mixed)
library(ggpubr)

#Representamos:
ggplot(finalSKI,aes(finalSKI$cg08884752 , finalSKI$Z_SCORE)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman")



#Sonda para el gen CALCB:cg16394018
#Necesito cargar la lista de Gencode que tiene el identificador del gen para hacer el genCALCB.

listaGenesGencode <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/listaGenesGencode.csv")

CosmicSS450K <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_result450KSS$Sample.Name,]

CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% listaGenesUSDC$nombreGenesUSDCSplitUnicos,]
#CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% listaGenesGencode$nombreGenesGencodeSplitUnlistUnicos,]

betas_result450KSSUnique <- betas_result450KSS[!duplicated(betas_result450KSS$Sample.Name),]

CosmicSS450KNombreGenes <- CosmicSS450KGenes$GENE_NAME
count <-0 
CosmicSS450KNombreGenes <- unique(CosmicSS450KNombreGenes)
for(i in CosmicSS450KNombreGenes){
  count <- count + 1
  assign(paste0("gen", i), CosmicSS450KGenes[CosmicSS450KGenes$GENE_NAME %in% i,])
}


expresionGen <- betas_result450KSSUnique[betas_result450KSSUnique$Sample.Name %in% genCALCB$SAMPLE_NAME ,c("Sample.Name", "cg16394018")]

finalCALCB <- merge(genCALCB, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")

plot <- ggplot(finalCALCB,aes(cg16394018,Z_SCORE)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

ggsave(plot, filename = "graficas/plotMetilacionvsExpresion.pdf", device = cairo_pdf, 
       width = 4, height = 3, units = "in")


DMPs0.15Signi$UCSC_RefGene_Name<- strsplit(as.character(DMPs0.15Signi$UCSC_RefGene_Name), ";")

DMPs0.15Signi2 <- unnest(DMPs0.15Signi,UCSC_RefGene_Name)

#DMPs0.15Signi2 <- unlist(DMPs0.15Signi,UCSC_RefGene_Name)


DMPs0.15Signi2 <- DMPs0.15Signi2[unique(DMPs0.15Signi2$UCSC_RefGene_Name),]

DMPs0.15Signi <- DMPs0.15Signi[!duplicated(DMPs0.15Signi$UCSC_RefGene_Name),]

DMPs0.15Signi2 <-    group_by(DMPs0.15Signi2$UCSC_RefGene_Name)



DMPs0.15Signi <- DMPs0.15Signi2


#automatizar:


sondasYnombres <- DMPs0.15Signi[,c("GencodeBasicV12_NAME","UCSC_RefGene_Name","Name","chr","pos")] #Añadir columan posicion y cromosoma
write.table(sondasYnombres, file="sondasYnombres.csv", sep=",", row.names=TRUE,  col.names=TRUE)


sondasYnombres <- sondasYnombres[!duplicated(sondasYnombres$GencodeBasicV12_NAME),]
sondasYnombres <- sondasYnombres[!duplicated(sondasYnombres$UCSC_RefGene_Name),]


#Estoy perdiendo 5 genes al coger solo el primero,debria de solucionarlo!!!!!!!!!!!!!!!!!!!!!!!!!!!
sondasYnombres <- DMPs0.15Signi[DMPs0.15Signi$UCSC_RefGene_Name %in%  listaUCSC,c("GencodeBasicV12_NAME","UCSC_RefGene_Name","Name")]
sondasYnombres <- unique(sondasYnombres$UCSC_RefGene_Name)

listaGenesGencode <- sapply(strsplit(DMPs0.15Signi$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])
#listaGenesGencode <- DMPs0.15Signi[,"GencodeBasicV12_NAME"]
listaGenesGencode <- listaGenesGencode[!duplicated(listaGenesGencode)]
listaGenesGencode <- na.omit(listaGenesGencode)

listaUCSC <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Name, ';', fixed=T), function(x) x[1])
listaUCSC <- listaUCSC[!duplicated(listaUCSC)]
listaUCSC <- na.omit(listaUCSC)

listaUCSC <- listaGenesUSDC$X


#Como no todos los genes de cosmic estan en mis betas results cojo los genes del betas resutls
nombreGenesBetasResult <- betas_result450K[5,6:ncol(betas_result450K)]



sondasYnombres$GencodeBasicV12_NAME <- sapply(strsplit(sondasYnombres$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])
sondasYnombres <- na.omit(sondasYnombres)

sondasYnombres$UCSC_RefGene_Name <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Name, ';', fixed=T), function(x) x[1])
sondasYnombres <- sondasYnombres[!duplicated(sondasYnombres$UCSC_RefGene_Name),]

DMPs0.15Signi$UCSC_RefGene_Name <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Name, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$UCSC_RefGene_Group <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Group, ';', fixed=T), function(x) x[1])

betas_result450KSS <- read.csv("betas_result_2021-08-06-13-40.csv", skip=7)

CosmicSSEPICSamples <- readRDS("CosmicSSEPICSamples.rds")



pdf("graficas expresion-metilacion/GraficasUCSC.pdf")

#Para ejecutar el for necesitamos cambiar los objetos: listaUCSC, CosmicSS450KGenes y sondasYnombres.
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

betas_result450KSS <- read.csv("betas_result_2021-08-06-13-40.csv", skip=7)

listaUCSC <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Name, ';', fixed=T), function(x) x[1])
listaUCSC <- listaUCSC[!duplicated(listaUCSC)]
listaUCSC <- na.omit(listaUCSC)

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")
CosmicSS450K <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_result450KSS$Sample.Name,]
CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% listaUCSC,]
i <- "ZNF608"
#Leemos una lista de  genes
for(i in listaUCSC){
  #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
  if(i %in% CosmicSS450K$GENE_NAME){
    
    j <- assign(paste0("gen", i), CosmicSS450K[CosmicSS450K$GENE_NAME %in% i,])
    
    #sonda <- sondasYnombres[sondasYnombres$UCSC_RefGene_Name %in% i, "Name" ]
    sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
    
    if(sonda  %in% colnames(betas_result450KSSUnique)){
      expresionGen <- betas_result450KSSUnique[betas_result450KSSUnique$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      
      
      plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman") + ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      print(plot)
    }
  }
}
dev.off()

ggsave(plot = lista , filename = "graficas/Jero.pdf" , device = cairo_pdf, 
       width = 4, height = 3, units = "in")

ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman") + geom_text(aes(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab,k$SAMPLE_NAME,"")),hjust = 1.1)


k2 <- k[,c("Z_SCORE",sonda)]
write.table(k2, file="ZNF608datosMetExp.csv", sep=",", col.names=TRUE,row.names=FALSE)



library(ggrepel)

ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman") + geom_text(label=k$SAMPLE_NAME,hjust = 1.1) 


ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman")  + geom_text_repel(label=k$SAMPLE_NAME) 

muestrasEnLab <- c("A-427","A-549","CAL-12T","Calu-1","Calu-3","EPLC-272H","H1264"
                   ,"H-1299","H1395","H-1437","H-1568","H157","H1623","H1703","H1975","H1993","H2106","H2170","H2172",
                   "H-2228","H-23","H-3122","H-3255","H-441","H-520","H-522","HCC-366","HCC-827","LCLC-103H","NCI-H2030","NCI-H460")


#Lineas celulares de pulmon:

'%!in%' <- function(x,y)!('%in%'(x,y))


cosmicMutaciones <- read.csv("V94_38_CLP_MUTANT.csv")

lineasCelulares <- CosmicSS450KGenes[CosmicSS450KGenes$SAMPLE_NAME %in% cosmicMutaciones$SAMPLE_NAME ,]
lineasNoCelulares <- CosmicSS450KGenes[CosmicSS450KGenes$SAMPLE_NAME %!in% cosmicMutaciones$SAMPLE_NAME ,]



#Preaparamos DMPsSingi para obtener las posiciones de los genes en las graficas.

DMPs0.15Signi$UCSC_RefGene_Name <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Name, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$UCSC_RefGene_Group <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Group, ';', fixed=T), function(x) x[1])

DMPs0.15Signi$GencodeBasicV12_NAME <- sapply(strsplit(DMPs0.15Signi$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_Group <- sapply(strsplit(DMPs0.15Signi$GencodeBasicV12_Group, ';', fixed=T), function(x) x[1])


#Heatmap con p-value
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

bVals <- readRDS("bVals.rds")

bVals0.15Signi <- bVals[match(rownames(DMPs0.15Signi),rownames(bVals)),]

bVals0.15Signi <- na.omit(bVals0.15Signi)

dim(bVals0.15Signi)
colnames(bVals0.15Signi) <- targets$EPIC.Sample.Name
heatmap3(bVals0.15Signi, ColSideCut=10)
heatmap3(bVals0.15Signi)

library(heatmap3)

targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")

heatmap3(bVals0.15Signi, ColSideCut=10, method = "ward.D2") 

heatmap3(bVals0.15Signi, ColSideCut=10, method = "single")  #"average", "mcquitty", "median" or "centroid".

heatmap3(bVals0.15Signi, ColSideCut=10, method = "complete") #distfun = function(x){dist(x, method = "euclidean")},     "complete", "average", "mcquitty", "median" or "centroid".

library(sjmisc) #str_contains
x <- c()
for (i in colnames(bVals0.15Signi)){
  if( str_contains(i, "P")){
    x <- c(x,"blue")
  }else{
    x <- c(x,"yellow")
  }
  
}

heatmap3(bVals0.15Signi, ColSideColors=x, ColSideLabs="Signature", Rowv = NULL,
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward.D", labCol= "")



BiocManager::install("pvclust")
library(pvclust)
pvclust(bVals0.15Signi)
heatmap3(pvclust(bVals0.15Signi))

#Ward2
data <- c(10,1,6,17)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.0006303



data
#Complete
data <- c(10,0,6,18)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 6.107e-05

#Heatmap unsupervised

BiocManager::install("sjmisc")
library(sjmisc)
library(gplots)
library(heatmap3)

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

#set seed

set.seed(1234)

bVals <- readRDS("bVals.rds")
bValsRSS <- bVals[sample(nrow(bVals), size=1000), ] #
saveRDS(bValsRSS,"bValsSS.rds")
bValsRSS <- readRDS("bValsSS.rds")

colnames(bValsRSS) <- targets$EPIC.Sample.Name

library(sjmisc) #str_contains
x <- c()
for (i in colnames(bValsRSS)){
  if( str_contains(i, "P")){
    x <- c(x,"blue")
  }else{
    x <- c(x,"yellow")
  }
  
}



heatmap3(bValsRSS, ColSideCut=1, ColSideColors=x, method = "ward.D2") 


#formato manel
#heatmap3(bValsRSS, ColSideColors=x, Rowv = NA, method = "ward.D2") 

#heatmap3(bValsRSS, ColSideColors=x, Rowv = NA, col=c("green","black", "red"), method = "ward.D2") 

#heatmap3(bValsRSS, ColSideColors=x, Rowv = NA,  col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward.D2", labCol= "") 

#heatmap3(bValsRSS, ColSideColors=x, Rowv = NA,scale = "row",  
# col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 


#buenas
#heatmap3(bValsRSS, ColSideColors=x, ColSideLabs="Signature",scale = "none", Rowv = NULL,
# col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward", labCol= "")


heatmap3(bValsRSS, ColSideColors=x, ColSideLabs="Signature", Rowv = NULL,
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward.D", labCol= "")



#heatmap.2(bValsRSS, ColSideColors=x, col = greenred(75), dendrogram="col",Rowv = NA,trace="none")

heatmap.2(bValsRSS, ColSideColors=x, col = greenred(75),trace="none", method = "ward.D")


head(rownames(bValsRSS))

#heatmap3(bValsRSS, ColSideColors=x, ColSideLabs="Signature", Rowv = NA,scale = "none",
col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward.D2", labCol= "",legendfun=function() showLegend(
  #legend = c("Group A", "Group B"),
  lwd = 3,
  cex = 1.1,
  #col = c("red", "blue"),
)) 

legend(x="right", legend=c("covid", "control"), fill=colorRampPalette(c("blue", "yellow"))(2))


heatmap3(bValsRSS, ColSideCut=1, ColSideColors=x, method = "centroid") 

#Ward2 1000 2º Ite
data <- c(11,5,5,13)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.0374

#Ward2 100
data <- c(11,6,5,12)
data <- matrix(data, ncol=2)
data
fisher.test(data)
p-value = 0.08441

#Ward2 100 2ª iteracion
data <- c(3,10,13,8)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.03865

#Complete 100 2ª iteracion
data <- c(4,10,12,8)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.09209

#Ward2 1000 3º Ite
data <- c(10,4,6,14)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.03488

#Ward2 1000 3º Ite
data <- c(12,5,4,13)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.03488




#Ward
data <- c(10,4,6,14)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.03488

data <- c(12,4,4,14)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.005



DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

bValsSS0.15 <- bVals[rownames(bVals) %in% rownames(DMPs0.15Signi),]

bVals0.15Signi <- bVals[match(rownames(DMPs0.15Signi),rownames(bVals)),]
colnames(bVals0.15Signi) <- targets$EPIC.Sample.Name

library(sjmisc) #str_contains
x <- c()
for (i in colnames(bVals0.15Signi)){
  if( str_contains(i, "P")){
    x <- c(x,"blue")
  }else{
    x <- c(x,"yellow")
  }
  
}


#heatmap3(bVals0.15Signi, ColSideCut=10, ColSideColors=x, method = "ward.D2") 


heatmap3(bVals0.15Signi, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 

#Ward
data <- c(11,1,5,17)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.0001806

heatmap3(bVals0.15Signi, ColSideColors=x, method = "ward.D", scale="none") 


#Ward
data <- c(10,0,6,18)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 6.107e-05

heatmap3(bVals0.15Signi, ColSideCut=10, ColSideColors=x, method = "complete") 

#Ward
data <- c(10,0,6,18)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 6.107e-05

#Ward
data <- c(10,0,6,18)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 6.107e-05


lista <- c()
count <- 0
#Pintar los genes join con cosmic(450K).
for(i in listaGenesGencode){
  if(i %in% CosmicSS450KNombreGenes){
    
    sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
    
    if(sonda  %in% colnames(betas_result450KSSUnique)){
      lista <- paste(lista, i)
      count<- count + 1
    }
  }
}
lista
dim(lista)



#Comprobar cuantas de mis sondas estan en 450K.
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sondas850K <- rownames(ann850k)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sondas450K <- rownames(ann450k)

sondas120 <- DMPs0.15Signi[,"Name"] #Añadir columan posicion y cromosoma

join <- sondas120[ sondas120 %in% sondas450K] #52 sondas

joinDMPS <- DMPs0.15Signi[ DMPs0.15Signi$Name %in% join,] #37 con genes

joinDMPS <- joinDMPS$UCSC_RefGene_Name
joinDMPS <- sapply(strsplit(joinDMPS, ';', fixed=T), function(x) x[1])

CosmicSS450KGenes$GENE_NAME

genesCosmic <- CosmicCompleteCL$GENE_NAME
genesCosmic <- unique(genesCosmic)

joinCosmic <- genesCosmic[ genesCosmic %in% joinDMPS]

#Unir UCSC y gencode
listaGenesGencode <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/listaGenesGencode.csv")
listaGenesUCSC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/listaGenesUSDC.csv")

listaGenesGencode <- listaGenesGencode$X
listaGenesUCSC <- listaGenesUCSC$X

joinGenes <- listaGenesUCSC[listaGenesUCSC %in% listaGenesGencode]
outerjoinGenes <- listaGenesUCSC[!listaGenesUCSC %in% listaGenesGencode]
outerjoinGenes2 <- listaGenesGencode[!listaGenesGencode %in% listaGenesUCSC]
joinGenes <- c(joinGenes,outerjoinGenes,outerjoinGenes2)
joinGenes

joinGenes<- as.data.frame(joinGenes)

joinGenes<-unique(joinGenes)
write.table(joinGenes, file="joinGenes.csv", sep=",", row.names=FALSE)

joinGeens2 <-  c(listaGenesUCSC,listaGenesGencode)
joinGeens2<-unique(joinGeens2)

joinGeens2[!joinGeens2 %in% joinGenes]

#Crear tabla con las 3 sondas
sondasInteres <- c("cg19028462","cg06055229","cg13493526","cg20033583")
betasInteres <- bVals[rownames(bVals) %in% sondasInteres,]
colnames(betasInteres) <- targets$EPIC.Sample.Name
write.table(betasInteres, file="betasInteres.csv", sep=",", row.names=TRUE)

library(heatmap3)
heatmap3(betasInteres, ColSideCut=10, method = "ward.D2") #, ColSideColors=x

tablaSondas <- DMPs0.15Signi[,c("Name","chr","pos","GencodeBasicV12_NAME","UCSC_RefGene_Name","P.Value")] #Añadir columan posicion y cromosoma
write.table(tablaSondas, file="tablaSondas.csv", sep=",", row.names=TRUE)


#Nos vamos a quedar con la anotacion de UCSC, y en las que no tenga nada le ponemos Gencode
sondasYnombres <- DMPs0.15Signi[,c("Name","GencodeBasicV12_NAME","UCSC_RefGene_Name")] #Añadir columan posicion y cromosoma
write.table(sondasYnombres, file="sondasYnombres.csv", sep=",", row.names=FALSE)


sondasYnombresUCSCFull <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/sondasYnombres.csv", TRUE, sep = ",")
sondasYnombresUCSCFull <- sondasYnombresUCSCFull$UCSC_RefGene_Name


sondasYnombresUCSCFullSplit<- strsplit(as.character(sondasYnombresUCSCFull), ";")

sondasYnombresUCSCFullSplitUnlist <-unlist(sondasYnombresUCSCFullSplit)

sondasYnombresUCSCFullSplitUnlistUnicos <- unique(sondasYnombresUCSCFullSplitUnlist)

nombreGenesGencode <- as.data.frame(nombreGenesGencodeSplitUnlistUnicos)

write.table(nombreGenesGencode, file="listaGenesGencode.csv", sep=",", row.names=TRUE)


#Graficas metilacion-expresion
setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

library(tidyverse)

BiocManager::install("sjmisc")
install.packages("dplyr")
library(rapport)
library(broom.mixed)
library(ggpubr)
library(sjmisc)


#Descargo los datos de metilacion de methylationDB:

betas_resultEPIC <- read.csv("betas_resultAllCellLinesEPIC373.csv", skip=7) #lung
#De las 120 sondas que introduzco solo encuentra 93.

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")


CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

CosmicSSEPIC <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name,]

joinGenes <- read.csv("joinGenes.csv")
joinGenes <- joinGenes$x


joinGenes2 <- joinGenes2$x

joinGenes2[!joinGenes2 %in% joinGenes]

CosmicSSEPICGenes <- CosmicSSEPIC[CosmicSSEPIC$GENE_NAME %in% joinGenes,]

sondasYnombres <- DMPs0.15Signi[,c("Name","GencodeBasicV12_NAME","UCSC_RefGene_Name")] #Añadir columan posicion y cromosoma
sondasYnombres$UCSC_RefGene_Name <- sapply(strsplit(sondasYnombres$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
sondasYnombres$GencodeBasicV12_NAME <- sapply(strsplit(sondasYnombres$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])


#CosmicSSEPICNombreGenes <- CosmicSSEPIC$GENE_NAME
#CosmicSSEPICNombreGenes <- unique(CosmicSSEPICNombreGenes)

#CosmicSSEPICSamples <- CosmicCompleteCL$SAMPLE_NAME[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name]
#CosmicSSEPICSamples <- unique(CosmicSSEPICSamples)
#saveRDS(CosmicSSEPICSamples, "CosmicSSEPICSamples.rds")

#for de EPIC.

i <- "MCF2L"

pdf("graficas expresion-metilacion/GraficasGencodeAllCellLinesEPIC373.pdf")

#Para ejecutar el for necesitamos cambiar los objetos: listaUCSC, CosmicSS450KGenes y sondasYnombres.
#Leemos una lista de  genes
for(i in joinGenes){
  #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
  if(i %in% CosmicSSEPICGenes$GENE_NAME){
    
    j <- assign(paste0("gen", i), CosmicSSEPICGenes[CosmicSSEPICGenes$GENE_NAME %in% i,])
    
    #sonda <- sondasYnombres[sondasYnombres$UCSC_RefGene_Name %in% i, "Name" ]
    sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
    
    sonda <- sonda[1]
    
    if(!is_empty(sonda) && sonda  %in% colnames(betas_resultEPIC)){
      expresionGen <- betas_resultEPIC[betas_resultEPIC$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      variablePromotor <- DMPs0.15Signi[i %in% DMPs0.15Signi$GencodeBasicV12_NAME,"GencodeBasicV12_Group"]
      
      
      
      plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman")+ geom_text(label=k$SAMPLE_NAME) + ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      print(plot)
    }
    #  lista<- append(lista,plot)
    #lista[count] <- plot                               # Store plots in list
    #count <- count + 1
    
  }
}
dev.off()

is.empty(sonda)
install.packages("sjmisc")

is_empty(sonda)
rm(x)


"cg05648629" %in% colnames(betas_resultEPIC)

ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman") + geom_text(label=k$SAMPLE_NAME)


#BLueprint

#1:Hay que ejecutarlo en linux.




#find.overlap Genomic Ranges:
setwd("/Users/jparra/Documents/ProyectoCovid/data/blueprint/")

DMPs0.15Signi <- read.csv("DMPs0.15Signi.csv")
#Anotamos DMPs para hg38



ss <- read.csv("ejemploSalidaBigWig.csv")

#DMPs0.15Signi$ranges <- paste(DMPs0.15Signi$pos,"-",DMPs0.15Signi$pos+2)

#Crear objeto GRanges
#GrangesObjectsCelulas <- GRanges(seqnames=ss$ï..Celulas, ranges=ss$Iranges, strand=ss$strand)


#GrangesObjectsMuestras <- GRanges(seqnames=rownames(DMPs0.15Signi), ranges=DMPs0.15Signi$pos, strand=DMPs0.15Signi$strand)





#Anotacion hg38
annhg38 <- read.csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7)

regiones38 <- annhg38[,c("Name","CHR_hg38","Start_hg38","End_hg38","Strand_hg38")]

# start_hg38 <- annhg38[annhg38$Name %in% DMPs0.15Signi$Name,c("Start_hg38","Name")]
#DMPs0.15SigniFinal <- merge(DMPs0.15Signi, start_hg38, by.x = "Name", by.y = "Name")

#End_hg38 <- annhg38[annhg38$Name %in% DMPs0.15Signi$Name,c("End_hg38","Name")]
#DMPs0.15SigniFinal <- merge(DMPs0.15SigniFinal, End_hg38, by.x = "Name", by.y = "Name")

#write.table(DMPs0.15SigniFinal, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/DMPs0.15SigniFinal.csv", sep=",", row.names=FALSE)
#DMPs0.15SigniFinal <- read.csv("bigwig/DMPs0.15SigniFinal.csv")

#GrangesObjectsCelulas <-makeGRangesFromDataFrame(ss, seqnames.field="chr", start.field="init", end.field="final", strand.field="strand")
#GrangesObjectsMuestras <- makeGRangesFromDataFrame(DMPs0.15SigniFinal, seqnames.field="chr", start.field="Start_hg38", end.field="End_hg38", strand.field="strand")

regiones38 <- na.omit(regiones38)

GrangesObjectsEPIC <- makeGRangesFromDataFrame(regiones38, seqnames.field="CHR_hg38", start.field="Start_hg38", end.field="End_hg38", strand.field="Strand_hg38")

# Save a single object to a file
saveRDS(GrangesObjectsEPIC, "GrangesObjectsEPIC.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DEL OBJETO RGSET EN RDS !!!!!!!!!!!!!
GrangesObjectsEPIC <- readRDS("GrangesObjectsEPIC.rds")
# Save a single object to a file
saveRDS( regiones38, " regiones38.rds")
# Restore it under a different name             !!!!!!!COMIENZO A PARTIR DEL OBJETO RGSET EN RDS !!!!!!!!!!!!!
regiones38 <- readRDS(" regiones38.rds")

bedGrapg <- read_bed_graph("/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/bedGraph/NeuFemale2.salida.bedGraph")

bedGrapgDF <- as.data.frame(bedGrapg)
head(bedGrapgDF)
find <- findOverlaps(GrangesObjectsEPIC,bedGrapg) 
find

#join <- sqldf("select * from bedGrapgDF f1 inner join DMPs f2 on (f1.start > f2.Start_hg38 and f1.start<= f2.End_hg38) ")



findDF<- as.data.frame(find)
findDF <-  findDF[!duplicated(findDF[,1]),]


sondas <-  regiones38[findDF[,1],"Name"]
sondas<- as.data.frame(sondas)
sondas$beta <- bedGrapgDF[findDF[,2],"score"]

write.table(sondas, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/NeuFemale2Full.csv", sep=",", row.names=FALSE)



monocito1 <- read.csv("bigwig/celulasEpiProfile/MonociteT0.csv")
monocito2 <- read.csv("bigwig/celulasEpiProfile/MonocyteNoneFemale.csv")
CD4 <- read.csv("bigwig/celulasEpiProfile/CD4Male.csv")
CD42 <- read.csv("bigwig/celulasEpiProfile/CD4Female329.csv")
CD43<- read.csv("bigwig/celulasEpiProfile/CD4Male5.2.csv")

CD4 <- read.csv("bigwig/celulasEpiProfile/CD4MaleFull.csv")
CD42 <- read.csv("bigwig/celulasEpiProfile/CD4Female329Full.csv")
CD43 <- read.csv("bigwig/celulasEpiProfile/CD4Male5.2Full.csv")
monocito1 <- read.csv("bigwig/celulasEpiProfile/MonociteT0Full.csv")
monocito2 <- read.csv("bigwig/celulasEpiProfile/MonocyteNoneMaleFull.csv")
monocito3 <- read.csv("bigwig/celulasEpiProfile/MonocyteNoneFemaleFull.csv")
CD8 <- read.csv("bigwig/celulasEpiProfile/CD8MaleFull.csv")
CD82 <- read.csv("bigwig/celulasEpiProfile/CD8FemaleFull.csv")
CD83 <- read.csv("bigwig/celulasEpiProfile/CD8MaleCBFull.csv")
NeuMale <- read.csv("bigwig/celulasEpiProfile/NeuMaleFull.csv")
NeuFemale <- read.csv("bigwig/celulasEpiProfile/NeuFemaleFull.csv")
NeuFemale2 <- read.csv("bigwig/celulasEpiProfile/NeuFemale2Full.csv")

celulas <- merge(CD4,CD42, by.x = "sondas", by.y = "sondas") #Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,CD43, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,monocito1, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,monocito2, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,monocito3, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,CD8, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,CD82, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,CD83, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,NeuMale, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,NeuFemale, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
celulas <- merge(celulas,NeuFemale2, by.x = "sondas", by.y = "sondas" ) #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge


colnames(celulas) <- c("sondas","CD4","CD42","CD43","monocito1","monocito2","monocito3","CD8","CD82","CD83","NeuMale","NeuFemale","NeuFemale2")
rownames(celulas) <- celulas$sondas
#rm(monocitos[,"sondas"])
celulas <- celulas[,2:13]
#heatmap3(monocitos)
#monocitos <-  monocitos[monocitos!=0,]
#monocitos <- monocitos[!is.na(monocitos),]
write.table(celulas, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/celulasYsondas.csv", sep=",", row.names=FALSE)
celulasYsondas <- read.csv("bigwig/celulasEpiProfile/celulasYsondas.csv")

colnames(celulas) <- c("CD4","CD42","CD43","monocito1")

#buscar valor por cpg:

CD4[CD4$sondas %in% "cg13493526",]
CD42[CD42$sondas %in% "cg13493526",]
CD43[CD43$sondas %in% "cg13493526",]
monocito1[monocito1$sondas %in% "cg13493526",]
monocito2[monocito2$sondas %in% "cg13493526",]
monocito3[monocito3$sondas %in% "cg13493526",]
CD8[CD8$sondas %in% "cg13493526",]
CD82[CD82$sondas %in% "cg13493526",]
CD83[CD83$sondas %in% "cg13493526",]
NeuMale[NeuMale$sondas %in% "cg13493526",]
NeuFemale[CD83$NeuFemale %in% "cg13493526",]
NeuFemale2[NeuFemale2$sondas %in% "cg13493526",]

# row_sub = apply(monocitos, 1, function(row) all(row !=0)) 
#monocitos <- monocitos[row_sub,] 
# row_sub = apply(monocitos, 1, function(row) all(row !=1)) 
# monocitos <- monocitos[row_sub,] 

#heatmap.2(as.matrix(celulas),scale="none",trace="none", col=greenred)
#Distancias y metodos de clustering
celulasSS <- celulas[sample(nrow(celulas), size=10000), ]

heatmap.2(as.matrix(celulasSS),scale="none",trace="none", col=greenred,dendrogram="column")
#heatmap3(as.matrix(celulasSS),scale="none",trace="none", col=greenred,dendrogram="column")


getwd()
bVals <- readRDS("bigWig/bVals.rds")
bValsDF <- as.data.frame(bVals)
bValsDF$sondas  <- rownames(bValsDF)
celulas$sondas <- rownames(celulas)
celulasYMuestras <- merge(celulas,bValsDF, by.x = "sondas", by.y = "sondas") #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
rownames(celulasYMuestras) <- celulasYMuestras$sondas
celulasYMuestras <- celulasYMuestras[,2:46]
write.table(celulasYMuestras, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/celulasYMuestras.csv", sep=",", row.names=TRUE)

celulasYMuestrasSS <- celulasYMuestras[sample(nrow(celulasYMuestras), size=10000), ]

heatmap.2(as.matrix(celulasYMuestrasSS),scale="none",trace="none", col=greenred,dendrogram="column")


#Covertura
bedGrapgCov <- read_bed_graph("/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/coverture/CD452Cov.salida.bedGraph")
bedGrapg <- read_bed_graph("/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/bedGraph/CD4Male5.2.salida.bedGraph")
bedGrapgDF <- as.data.frame(bedGrapg)
head(bedGrapgDF)
bedGrapgCovDF <- as.data.frame(bedGrapgCov)
head(bedGrapgCovDF)
find <- findOverlaps(bedGrapgCov,bedGrapg) 
find
lecturaYcovertura <- merge(bedGrapgDF,bedGrapgCovDF, by.x = "start", by.y = "start") 
#colnames(lecturaYcovertura)
lecturaYcovertura <- lecturaYcovertura[,c("seqnames.x","start","end.x","score.x","score.y")]
lecturaYcoverturaCov <- subset(lecturaYcovertura,score.y>10)


GrangeslecturaYcovertura <- makeGRangesFromDataFrame(lecturaYcoverturaCov2, seqnames.field="seqnames.x", start.field="start", end.field="end.x", strand.field="strand.x")

find <- findOverlaps(GrangesObjectsEPIC,GrangeslecturaYcovertura) 

find
findDF<- as.data.frame(find)
findDF <-  findDF[!duplicated(findDF[,1]),]

nrow(subset(lecturaYcovertura, score.x!=1))
nrow(subset(lecturaYcovertura, score.x!=1 | score.x!=0))


lecturaYcoverturaCov2 <- subset(lecturaYcovertura, score.x>=0.1 & score.x<=0.9)

write.table(lecturaYcoverturaCov, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/lecturaYcoverturaCD452.csv", sep=",", row.names=TRUE)

lecturaYcovertura <- read.csv("bigwig/celulasEpiProfile/lecturaYcoverturaMonociteT0.csv")


sondas <-  regiones38[findDF[,1],"Name"]
sondas<- as.data.frame(sondas)
sondas$beta <- lecturaYcoverturaCov2[findDF[,2],"score.x"]

write.table(sondas, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/MonociteT0Cov.csv", sep=",", row.names=FALSE)


MonociteT0 <- read.csv("bigwig/celulasEpiProfile/MonociteT0Cov.csv")
CD452 <- read.csv("bigwig/celulasEpiProfile/CD452Cov.csv")

celulas <- merge(MonociteT0,CD452, by.x = "sondas", by.y = "sondas")
celulasYMuestras <- merge(sondas,bValsDF, by.x = "sondas", by.y = "sondas") #   ,all=TRUE    Quiza merged que ponga 0 en las que no hacen  merge
rownames(celulasYMuestras) <- celulasYMuestras$sondas
celulasYMuestras <- celulasYMuestras[,2:36]

celulasYMuestrasSS <- celulasYMuestras[sample(nrow(celulasYMuestras), size=10000), ]

heatmap.2(as.matrix(celulasYMuestrasSS),scale="none",trace="none", col=greenred,dendrogram="column")
write.table(celulasYMuestras, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/celulasYMuestras0.1.csv", sep=",", row.names=FALSE)
celulasYMuestras <- read.csv("bigwig/celulasEpiProfile/celulasYMuestras0.1.csv")


#Eliminar todos los valores 1 o 0
lecturaYcoverturaCov2 <- subset(lecturaYcovertura, score.x>=0.1 & score.x<=0.9)



#Episcore


remotes::install_github("immunogenomics/presto")
BiocManager::install("org.Hs.eg.db")
devtools::install_github("aet21/EpiSCORE")
library(EpiSCORE)
library(presto)
BiocManager::install("EpiDISH")
library(EpiDISH)

setwd("/Users/jparra/Documents/ProyectoCovid/data/blueprint/")

help(EpiSCORE)
??EpiSCORE 

#Estudiar como hace el match!!!!!!!! por nombre de gen??

## Building the expression reference matrix

data(lungSS2mca1)
#Contains the normalized SmartSeq2 lung atlas from MCA1 with rows
#labeling 15845 genes and columns labeling 1641 cells, an integer cell-type index vector, and the names of the cell-types
#Tenemos 3 objetos, una matriz con los valores de expresion para 1641 celulas de 4 tipos diferentes. Y dos vectores con lo tipos de celulas y el indice
ncpct.v <- summary(factor(celltypeSS2.idx));
names(ncpct.v) <- celltypeSS2.v;
print(ncpct.v);

# data(lung10Xmca1) Validation


expref.o <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH=rep(3,4))
print(dim(expref.o$ref$med))
head(expref.o$ref$med)

estF.m <- epidish(lung10Xmca1.m,ref.m=expref.o$ref$med,method="RPC",maxit=200)$estF

par(mfrow=c(1,4));
par(mar=c(4,4,2,1));
for(ct in 1:4){
  boxplot(estF.m[,ct] ~ celltypeSS2.v[celltype10X.idx],ylab="EstFrac",xlab="10X cell-type",main=paste("EstFrac-",celltypeSS2.v[ct],sep=""),pch=23,ylim=c(0,1));
}

pred.idx <- apply(estF.m,1,which.max);
acc <- length(which(pred.idx==celltype10X.idx))/length(celltype10X.idx);
print(paste("Overall accuracy=",round(acc,2),sep=""));

## Building the DNAm reference matrix


refMscm2.m <- ImputeDNAmRef(expref.o$ref$med,db="SCM2",geneID="SYMBOL");
refMrmap.m <- ImputeDNAmRef(expref.o$ref$med,db="RMAP",geneID="SYMBOL");

refMmg.m <- ConstMergedDNAmRef(refMscm2.m,refMrmap.m)
print(dim(refMmg.m));
head(refMmg.m);

load("dataExampleLung.rda")
print(dim(avSIM.m));
#avSIM son los betaValues de 100 mezclas de celulas con el gen entrez id del promotor.
#avSIM.m: 100 simulated DNAm mixtures with DNAm values summarized at the gene promoter level.
#trueW.m: true mixing fractions for the 100 mixtures.
# avLUSCtss.m: the TCGA LUSC DNAm dataset of 316 samples, summarized at the gene-promoter level.
#phenoLUSC.lv: list of phenotype values for the TCGA dataset.


plot(density(refMmg.m[,5]),lwd=2,xlab="Weight",main="");
abline(v=0.4,lwd=2,col="red");

print(paste("Number of selected genes=",length(which(refMmg.m[,5]>0.4)),sep=""));
estF.o <- wRPC(data=avSIM.m,ref=refMmg.m,useW=TRUE,wth=0.4,maxit=200);

print(cor(estF.o$estF,trueW.m));
pcc.v <- diag(cor(estF.o$estF,trueW.m));

par(mfrow=c(1,4));
par(mar=c(4,4,2,1));
for(ct in 1:4){
  plot(estF.o$estF[,ct],trueW.m[,ct],pch=23,xlim=c(0,1),ylim=c(0,1),main=colnames(trueW.m)[ct],xlab="fCT(Estimated)",ylab="fCT(True)",cex=0.5);
  abline(a=0,b=1,col="green",lty=2,lwd=2);
  text(x=.5,y=0.9,paste("PCC=",round(pcc.v[ct],2),sep=""),font=2,cex=1); 
}

#Real data:
estF.o <- wRPC(avLUSCtss.m,ref=refMmg.m,useW=TRUE,wth=0.4,maxit=200);
pv <- wilcox.test(estF.o$est[,"Epi"] ~ phenoLUSC.lv$Cancer,alt="less")$p.value;
print(pv)


# Having obtained the cell-type fractions for the main cell-types in all lung-tissue samples, 
# we can now aim to identify differentially methylated cell-types (DMCTs) in cancer,
# i.e. cell-type specific differentially methylated cytosines. Since we have a DNAm reference matrix,
# we can apply a reference-based method like CellDMC @Zheng2018,
# which incorporates statistical interaction terms between phenotype (normal vs cancer) 
# and estimated cell-type fractions. Because the input are the estimated cell-type fractions themselves,
# CellDMC is run **over all CpGs**, that is,
# while the inference of cell-type fractions was done by consider promoter DNAm levels,
# identifying DMCTs is done at single cytosine resolution level.
# `CellDMC` is part of the `EpiDISH` package loaded in earlier and we provide detailed vignette in 
# that package to explain how it is run. Because the full LUSC dataset `bmiqLSCCrmRS.m` is too large for 
# inclusion here, we just display the syntax:




cdmc.o <- CellDMC(bmiqLSCCrmRS.m,pheno.v=phenoLUSC.lv$Cancer+1, frac.m=estF.o$est, adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, sort = FALSE, mc.cores = 4);

dmctLUSC.lv <- list();
for(ct in 1:4){
  dmctLUSC.lv[[ct]] <- rownames(cdmc.o$coe[[ct]][which(cdmc.o$dmct[,1+ct]!=0),]);
}




BiocManager::install("SuperExactTest")
library(SuperExactTest);
res.o <- supertest(dmctLUSC.lv,n=395963)
plot(res.o,"landscape",sort.by="size")



#Jero
setwd("/Users/jparra/Documents/ProyectoCovid/data/blueprint/")

bVals <- readRDS("bigWig/bVals.rds")

celulas <- read.csv("bigwig/celulasEpiProfile/celulas")


#cdmc.o <- CellDMC(bVals,pheno.v=bVals$205061410147_R03C01 , frac.m=estF.o$est, adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, sort = FALSE, mc.cores = 4);

estF.o <- wRPC(data=bVals,ref=refMmg.m,useW=FALSE);

celulas$weight <- celulas$beta

estF.o <- wRPC(data=bVals,ref=celulas)
#estF.o <- wRPC(avLUSCtss.m,ref=refMmg.m,useW=TRUE,wth=0.4,maxit=200);


rownames(bVals) <- (1:865859)
rownames(celulas) <- (1:63)

bValsSS <- bVals[1:50,1:5]
phenoJero <- c("0","1","1","1","1","1","1","1","1","1","1","1","0","0","0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1")


write.table(bValsSS, file="/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/bValsSSPEC.csv", sep=",", col.names=TRUE)
setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")


bVals <- as.matrix(bVals)
celulas <- as.matrix(celulas)

bValsSS <- as.matrix(bValsSS)
colnames(bVals) <- targets$EPIC.Sample.Name

celldmc.o <- CellDMC(bVals, phenoJero, frac.m) 

out.l <- epidish(DummyBeta.m, celulas, method = 'RPC')
out.l <- epidish(bVals, celulas, method = 'RPC')

out.l <- epidish(bValsSS, celulas, method = 'RPC')

#Examples:

data(centEpiFibIC.m)
data(DummyBeta.m)
out.l <- epidish(DummyBeta.m, centEpiFibIC.m, method = 'RPC')
frac.m <- out.l$estF
pheno.v <- rep(c(0, 1), each = 5)
celldmc.o <- CellDMC(DummyBeta.m, pheno.v, frac.m) 
celldmc.o$dmct

###########

lung <-  load("C:/Users/jparra/Downloads/LungRef.rda")
lung <-  load("/Users/jparra/Downloads/dataExampleLung.rda")

lung <- data(dataExampleLung)

class(lung)
head(lung)
str(lung)
lung@
  head(celldmc.o$coe)
class(centEpiFibIC.m)
is.numeric(celulas)
rownames() %in% rownames(centEpiFibIC.m)
rownames(centEpiFibIC.m) %in% rownames(DummyBeta.m)
table( rownames(centEpiFibIC.m) %in% rownames(DummyBeta.m))


write.table(celulas, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/celulasSS.csv", sep=",", row.names=FALSE)

celulas <- read.csv("bigWig/celulasEpiProfile/celulas.csv")


bValsSS <- unique(bValsSS)
celulas <- unique(celulas)
class(celulas)
unique.pairs <- unique( bValsSS[ duplicated(bValsSS) ,])
grps <- apply( bValsSS, 1, function(x) if ( any( duplicated( rbind(unique.pairs, x))) ) { paste(x[1],x[2], sep="->")}else{NA} )
grps

out.l <- epidish(DummyBeta.m, celulas, method = 'RPC')
out.l <- epidish(bValsSS, centEpiFibIC.m, method = 'RPC')

BiocManager::install("paste.matrix")
library(paste.matrix);

bind <- rbind(bValsSS,celulas)
bind <- paste.matrix(bValsSS,celulas)



cdmc.o <- CellDMC(bmiqLSCCrmRS.m,pheno.v=phenoLUSC.lv$Cancer+1, frac.m=estF.o$est, adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, sort = FALSE, mc.cores = 4);




install.packages('MASS')
library(MASS)
#s in seq_len(ncol(DummyBeta.m))

if (qr(DummyBeta.m)$rank < ncol(DummyBeta.m)) 
  
  
  rlm.o <- rlm(bValsSS[, 5] ~ celulas)
DummyBeta.m[, 1]
dim(celulas)

map.idx <- match(rownames(bVals), rownames(celulas))

beta.m <- bVals
ref.m <- celulas
maxit <- 20
s <- 1
#RPC: funcion de epidish
DoRPC <- function(beta.m, ref.m, maxit) {
  map.idx <- match(rownames(ref.m), rownames(beta.m))
  rep.idx <- which(is.na(map.idx) == FALSE)
  data2.m <- beta.m[map.idx[rep.idx], ]
  ref2.m <- ref.m[rep.idx, ]
  est.m <- matrix(nrow = ncol(data2.m), ncol = ncol(ref2.m))
  colnames(est.m) <- colnames(ref2.m)
  rownames(est.m) <- colnames(data2.m)
  for (s in seq_len(ncol(data2.m))) {
    rlm.o <- rlm(data2.m[, s] ~ ref2.m, maxit = maxit)
    coef.v <- summary(rlm.o)$coef[2:(ncol(ref2.m) + 1), 1]
    coef.v[which(coef.v < 0)] <- 0
    total <- sum(coef.v)
    coef.v <- coef.v/total
    est.m[s, ] <- coef.v
  }
  return(list(estF = est.m, ref = ref2.m, dataREF = data2.m))
}

out.o <- DoRPC(bValsSS, celulas,20)

out.o <- DoRPC(DummyBeta.m, centEpiFibIC.m,20)

qr(bValsSS)$rank < ncol(bValsSS)

out.o<- list(estF = est.m, ref = ref2.m, dataREF = data2.m)
out.o$estF

rlm.o <- funcion(bValsSS[, 1] ~ celulas, maxit = 20)
rlm.o <- rlm(bValsSS[, 1] ~ celulas, maxit = 20)

write.table(out.o$estF, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/deconvolution.csv", sep=",", col.names=TRUE)


targets <- read.metharray.sheet("/Users/jparra/Documents/ProyectoCovid/data/blueprint", pattern="SampleSheetAutop34.csv")

bVals <- readRDS("bigWig/bVals.rds")
colnames(bVals) <- targets$EPIC.Sample.Name
avDNAm.m <- constAvBetaTSS(bVals,type="850k")

estF.o <- wRPC(avDNAm.m,ref=refMmg.m,useW=TRUE,wth=0.2,maxit=200);
estF.o$estF
estF.o$ref

write.table(estF.o$estF, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/deconvolution2.csv", sep=",", col.names=TRUE)


#pv <- wilcox.test(estF.o$estF[,"IC"] ~ phenoLUSC.lv$Cancer,alt="less")$p.value;
pv <- wilcox.test(estF.o$estF[,"IC"] ~ targets$TISSUE,alt="less")

print(pv)
#La concentracion de celulas del sistema inmune es estadisticamente mas abundante en los pacientes con covid.



cdmc.o <- CellDMC(avDNAm.m,pheno.v=targets$TISSUE, frac.m=estF.o$estF, adjPMethod = "fdr", adjPThresh = 0.6, cov.mod = NULL, sort = FALSE, mc.cores = 1);
cdmc.o$coe$IC
cdmc.o$dmct

cdmc.o2 <- CellDMC(bVals,pheno.v=targets$TISSUE, frac.m=estF.o$estF, adjPMethod = "fdr", adjPThresh = 0.6, cov.mod = NULL, sort = FALSE, mc.cores = 1);

ct <- 1
dmctLUSC2.lv <- list();
for(ct in 1:4){
  dmctLUSC2.lv[[ct]] <- rownames(cdmc.o2$coe[[ct]][which(cdmc.o2$dmct[,1+ct]!=0),])
}

dmctLUSC.lv <- list();
for(ct in 1:4){
  dmctLUSC.lv[[ct]] <- rownames(cdmc.o$coe[[ct]][which(cdmc.o$dmct[,1+ct]!=0),])
}
cdmc.o2$dmct[,5]!=0

dmctLUSC.lv[4]

BiocManager::install("SuperExactTest")
library(SuperExactTest);

res.o <- supertest(dmctLUSC.lv,n=395963);
plot(res.o,"landscape",sort.by="size");





ggplot(a, aes(x = cellType, y = value, fill = group)) +    # Create boxplot chart in ggplot2
  geom_boxplot() + theme_bw()










funcionRLM <- function ....


update.packages()


#MAtrices de referencia del articulo episcore:
setwd("/Users/jparra/Documents/ProyectoCovid/data/blueprint/")

Fibro <- read_rds("Alveolar.Cellview.Rds")
Fibro <- read_rds("FibroCellview.Rds")

Fibro <- readRDS("GrangesObjectsEPIC.rds")
unzip("E-MTAB-6149.processed.1.zip")
Fibro <- readRDS("Fibro.Cellview.Rds")


temp <- readRDS.gz('Alveolar.Cellview.Rds')
temp <- readRDS.gz('FibroCellview.rds')

Fibro <- read_rds("GrangesObjectsEPIC.ds")

load("E-MTAB-6149.processed.1")

load("Alveolar.Cellview.Rds")




#Human Cell Atlas

#install.packages("Seurat")
#install.packages("spatstat.core")
library(spatstat.core)
library(Seurat)
library(readr)


rm(seurat)

seurat <- read_rds("/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/COVID-19 severity correlates with airway epithelium/covid_nbt_main.rds")
seurat@
plot1 <- VariableFeaturePlot(Fibro)


seurat[["RNA"]]@scale.data
all.genes <- rownames(seurat)
Fibro <- ScaleData(seurat, features = all.genes)
#Error: cannot allocate vector of size 27.2 Gb

seuratData <- ScaleData(seurat)




#Cell Sequence and Cell Label (pmbc is my data)

write.table(seurat@active.ident, file='Convert_UMI_Label.tsv', quote=FALSE, sep='\t', col.names = TRUE)

ident <- seurat@active.ident
identDF <- as.data.frame(ident)
identUnique <- unique(identDF$ident)

#Gene counts per cell

write.table(pbmc@assays[["RNA"]]@counts, file='Gene_Count_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)

view(seurat@assays[["RNA"]]@counts)


counts <- seurat@assays[["RNA"]]@counts
counts@Dimnames[1]

misGenes <- counts@Dimnames[1] %in% joinGenes
misGenes <- joinGenes  %in% counts@Dimnames[1]


genesSeurat <- counts@Dimnames[1]

genesSeurat <- unlist(counts@Dimnames[1])

misGenes <- joinGenes  %in% genesSeurat

#misGenes <- "RABGGTB"  %in% genesSeurat

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

joinGenes <- read.csv("joinGenes.csv")
joinGenes <- joinGenes$x

Fibro <- ScaleData(seurat, features = joinGenes)


# How can I extract expression matrix for all NK cells (perhaps, to load into another package)
nk.raw.data <- as.matrix(GetAssayData(seurat, slot = "counts")[, WhichCells(seurat, ident = "NK")])  #64 celulas

raw.data <- as.matrix(GetAssayData(seurat, slot = "counts")[,])

nkNeu.raw.data <- as.matrix(GetAssayData(seurat, slot = "counts")[, WhichCells(seurat, ident = c("NK","Neu"))])  #Neu: 32553



# How many cells are in each cluster
table(Idents(seurat))

str(Idents(seurat))

DimPlot(seurat, reduction = "umap")

DimPlot(seurat, reduction = "umap",raster=FALSE)



seurat$CellType <- Idents(seurat)
seurat2 <- seurat


Idents(seurat2) <- "severity"

seurat2$severity


DimPlot(seurat2, reduction = "umap")

seurat2$CellType <- Idents(seurat)


Idents(seurat) <- "CellType"

prop.table(table(Idents(seurat2)))

Idents(pbmc)





cluster.averages2 <- AverageExpression(seurat)
head(cluster.averages2[["RNA"]][, 1:24])

expresion <- cluster.averages[["RNA"]][, 1:24]

expresion <- as.data.frame(expresion)
expresionMisGenes <- expresion[rownames(expresion) %in%  joinGenes,]


library(gplots) #heatmap.2

#heatmap.2(as.matrix(expresionMisGenes),scale="none",trace="none", col=greenred,dendrogram="column")

#Depmap5$gen <- zscore(Depmap5$gen, dist="gamma", shape=0.5)

library(limma)

#expresionMisGenesZscore <- zscore(expresionMisGenes$Basal , dist="gamma", shape=0.5)
#expresionMisGenesZscore <- as.data.frame(expresionMisGenesZscore)
#expresionMisGenesZscore$`Secretory-diff` <- zscore(expresionMisGenes$`Secretory-diff` , dist="gamma", shape=0.5)
#expresionMisGenesZscore$Secretory <- zscore(expresionMisGenes$Secretory , dist="gamma", shape=0.5)


expresionMisGenesZscore <- zscore(as.matrix(expresionMisGenes) , dist="gamma", shape=0.5)


expresionMisGenesZscore[expresionMisGenesZscore == '-Inf'] <- 0

heatmap.2(expresionMisGenesZscore,scale="none",trace="none", col=greenred,dendrogram="column")



table(Idents(seurat), seurat$severity)


#Idents(seurat) <- "severity"


#WhichCells(seurat, idents = "control")

# How can I extract expression matrix for all control cells (perhaps, to load into another package)
#control.raw.data <- as.matrix(GetAssayData(seurat, slot = "counts")[, WhichCells(seurat, ident = "control")])


#expresionMisGenesControl <- control.raw.data[rownames(control.raw.data) %in%  joinGenes,]


library(readr)

seurat <- read_rds("/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/COVID-19 severity correlates with airway epithelium/covid_nbt_main.rds")

seuratControl <- subset(seurat, subset = severity == "control")

Idents(seuratControl)
Idents(seurat)

table(Idents(seurat), seurat$severity)
table(Idents(seuratControl), seuratControl$severity)

seuratCovid <- subset(seurat, subset = severity == "critical")
table(Idents(seuratCovid), seuratCovid$severity)

seuratCovid


Idents(seuratCovid)
Idents(seuratControl)



#IRF2BPL


cluster.averagesControl <- AverageExpression(seuratControl)
head(cluster.averagesControl[["RNA"]][, 1:19])

seuratControlExpresion <- cluster.averagesControl[["RNA"]][, 1:19]

seuratControlExpresion <- as.data.frame(seuratControlExpresion)
seuratControlExpresionIRF2BPL <- seuratControlExpresion[rownames(seuratControlExpresion) %in%  "IRF2BPL",]

seuratControlExpresionIRF2BPLZscore <- zscore(as.matrix(seuratControlExpresionIRF2BPL) , dist="gamma", shape=0.5)

seuratControlExpresionIRF2BPLZscore <- as.data.frame(seuratControlExpresionIRF2BPLZscore)


rownames(seuratControlExpresionIRF2BPLZscore) <- "CN"





cluster.averagesCovid <- AverageExpression(seuratCovid)
head(cluster.averagesCovid[["RNA"]][, 1:24])


seuratCovidExpresion <- cluster.averagesCovid[["RNA"]][, 1:24]


seuratCovidExpresion <- as.data.frame(seuratCovidExpresion)
seuratCovidExpresionIRF2BPL <- seuratCovidExpresion[rownames(seuratCovidExpresion) %in%  "IRF2BPL",]

seuratCovidExpresionIRF2BPLZscore <- zscore(as.matrix(seuratCovidExpresionIRF2BPL) , dist="gamma", shape=0.5)

seuratCovidExpresionIRF2BPLZscore <- as.data.frame(seuratCovidExpresionIRF2BPLZscore)


rownames(seuratCovidExpresionIRF2BPLZscore) <- "CV"





seuratControlExpresionIRF2BPLZscore[seuratControlExpresionIRF2BPLZscore == "-Inf"] <- 0

#seuratControlExpresionIRF2BPLZscore <- seuratControlExpresionIRF2BPLZscore[seuratControlExpresionIRF2BPLZscore != 0] 


colnames(seuratControlExpresionIRF2BPLZscore)

seuratCovidExpresionIRF2BPLZscore <- seuratCovidExpresionIRF2BPLZscore[colnames(seuratCovidExpresionIRF2BPLZscore) %in% colnames(seuratControlExpresionIRF2BPLZscore)] 

newdf <- rbind(seuratControlExpresionIRF2BPLZscore, seuratCovidExpresionIRF2BPLZscore)

library(gplots)
heatmap.2(as.matrix(newdf),scale="none",trace="none", col=greenred,dendrogram="column")

#La expresion del gen IRF2BPL en citotoxic cells (CTL) es bastante mayor en covid  que los en controles
#No hay ningun neutrofilo en los controles


boxplot(rownames(seuratControlExpresionIRF2BPLZscore) ~ colnames(seuratControlExpresionIRF2BPLZscore), seuratCovidExpresionIRF2BPLZscore)

boxplot(newdf$IFNResp,)


library(ggplot2)

ggplot(newdf, aes(x = Basal , y = CV , fill = Basal) + geom_boxplot())


ggplot(newdf, aes_(x = newdf$IFNResp , y = newdf$`Secretory-diff` , fill = newdf$IFNResp) + geom_boxplot())






identificadoresSeurat <- colnames(control.raw.data)





#COVID LUNG. seurat a partir de gen expresion

genExpresion <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/GSE171524_processed_data.csv")

joinGenes <- read.csv("joinGenes.csv")
genExpresionSS <- genExpresion[genExpresion$X %in% joinGenes,1:2000 ]

genExpresionSS <- genExpresion[,1:2000 ]


rownames(genExpresionSS) <- genExpresionSS$X
genExpresionSS <- genExpresionSS[,2:2000]


saveRDS(genExpresionSS, "/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/genExpresionSS.rds")
genExpresion <- readRDS("/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/genExpresion.rds")

saveRDS(genExpresionSS, "/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/genExpresionSS.rds")
genExpresionSS <- readRDS("/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/genExpresionSS.rds")

identificadores <- colnames(genExpresion)


identificadores %in% identificadoresSeurat

length(identificadores)
identificadores <- identificadores[2:length(identificadores)]

identificadores <- substr(identificadores,1,16)


identificadoresSeurat <- substr(identificadoresSeurat,1,16)

identificadores %in% identificadoresSeurat

identificadoresEnComun <- identificadores[identificadores %in% identificadoresSeurat]

FindAllMarkers(genExpresionSS)

FindAllMarkers(seurat)









#create seurat
install.packages('SeuratObject')
library(SeuratObject)
seurat <- CreateSeuratObject(counts = genExpresionSS, project = "jero", min.cells = 1, min.features = 1)

#seurat <- CreateSeuratObject(counts = genExpresionSS, project = "jero", names.field = 3, min.cells = 1, min.features = 1)


FindAllMarkers(seurat)

seurat2 <- CreateSeuratObject(counts = genExpresion, project = "jero", min.cells = 1, min.features = 2)

FindAllMarkers(seuratBueno)

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)








#pipeline: https://satijalab.org/seurat/archive/v3.1/pbmc3k_tutorial.html



GSM <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/GSM5226592_L12cov_raw_counts.csv")

rownames(GSM) <- GSM$X

GSM <- GSM[,2:length(GSM)]


seurat <- CreateSeuratObject(counts = GSM, project = "jero", min.cells = 1, min.features = 2)


seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")


VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(seurat), 10)

all.genes <- rownames(seurat)

seurat <- ScaleData(seurat, features = all.genes)


seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))


DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(seurat, dims = 1:4, cells = 500, balanced = TRUE)


seurat <- JackStraw(seurat, num.replicate = 100)

seurat <- ScoreJackStraw(seurat, dims = 1:25)

JackStrawPlot(seurat, dims = 1:30)


ElbowPlot(seurat)


saveRDS(seurat, file = "/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/seurat.rds")
seurat<- readRDS("/Users/jparra/Documents/ProyectoCovid/data/HumanCellAtlas/molecular single-cell lung atlas/seurat.rds")

library(Seurat)
seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 0.5)


head(Idents(seurat), 5)

seurat <- RunUMAP(seurat, dims = 1:20)

DimPlot(seurat, reduction = "umap")



FeaturePlot(seurat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                 "CD8A"))


FeaturePlot(seurat, features = "IRF2BPL")

FeaturePlot(seurat, features = "CALCB")



#Cibersort

#Vamos a construir la matriz de las muestras para cibersort con una muestra por fila y en las columnas los genes,
setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")
dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")

bVals<- readRDS("bVals.rds")
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

bValsSS <- bVals[rownames(bVals) %in% DMPs0.15Signi$Name,]
colnames(bValsSS) <- targets$EPIC.Sample.Name

bValsT <- t(bValsSS)


bValsDF <- as.data.frame(bVals)
#colnames(bValsSS) <- targets$EPIC.Sample.Name
bValsT <- t(bValsSS)
#bValsDFT <- as.data.frame(bValsT)

ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k$UCSC_RefGene_Name <- sapply(strsplit(ann850k$UCSC_RefGene_Name, ';', fixed=T), function(x) x[1])

colnames(bValsSS) <- ann850k$UCSC_RefGene_Name[ann850k$Name %in% colnames(bValsSS)]
rownames(bValsSS) <- ann850k$UCSC_RefGene_Name[ann850k$Name %in% rownames(bValsSS)]

#bValsDFT <- as.data.frame(bValsT)

bValsSSDF <- as.data.frame(bValsSS)

bValsT <- bValsSSDF[ -which(rownames(bValsT) %in% "") ,]

bValsT <- bValsSSDF[ -which(str_contains(rownames(bValsSSDF), "NA.")) ,]


BiocManager::install("sjmisc")
library(sjmisc)
str_contains(rownames(bValsSSDF), "NA.")




bValsT <- bValsT[ , -which(colnames(bValsT) %in% NA)]
bValsSS <- t(bValsT)


write.table(bValsSS, file="bValsSS.csv", sep=",", row.names=TRUE, col.names=TRUE)


bValsDF$genes <- ann850k$UCSC_RefGene_Name[ann850k$Name %in% rownames(bVals)]

bValsDF$genes <- sapply(strsplit(bValsDF$genes, ';', fixed=T), function(x) x[1])

rownames(bValsDF) <- bValsDF$genes
rm(bValsDF$genes)

bValsDF <- t(bValsDF)


bValsSSDF <- as.data.frame(bValsSS)
head(bValsDF)







devtools::install_github("IOBR/IOBR")
library(IOBR)
#BiocManager::install("ComplexHeatmap") version incompatible, marcar en la ventana de paquetes antes de devtools


cibersortAbsolutNoSum<-CIBERSORT(sig_matrix = lm22, mixture_file = bValsSSDF, perm = 1000, QN=TRUE, absolute=TRUE, abs_method = "no.sumto1")
write.table(cibersortAbsolutNoSum, file="cibersortAbsolutNoSum.csv", sep=",", row.names=TRUE, col.names=TRUE)



#BiocManager::install("MethylCIBERSORT")
#install.packages(c("caret", "glmnet", "NMF"))



#install.packages("MethylCIBERSORT")
devtools::install_local("MethylCIBERSORT")
library(MethylCIBERSORT)
??MethylCIBERSORT
data("V2_Signatures")
knitr::kable(head(Signatures$haematopoietic_neoplasm.other_v2_Signature.txt))

Signatures$haematopoietic_neoplasm.other_v2_Signature.txt




load("./DiscoveryCohort.rda")
data("StromalMatrix_V2")


Int <- intersect(rownames(bValsDF), rownames(Stromal_v2))
Mat <- bValsDF[match(Int, rownames(bValsDF)),]
Stromal_v2 <- Stromal_v2[match(Int, rownames(Stromal_v2)),]

RefData <- Stromal_v2
RefPheno <- Stromal_v2.pheno


Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = FALSE,
                              export = TRUE,
                              sigName = "MyReference",
                              Stroma.matrix = RefData,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = RefPheno)


Prep.CancerType(Beta = Mat, Probes = rownames(Signature$SignatureMatrix), fname = "MixtureMatrix")





#Hacer metilacion-expresion con lineas celulares sanguineas.

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

library(tidyverse) #read_tsv
library(ggrepel)
library(data.table)

#library(rapport)
#library(broom.mixed)
library(ggpubr)
library(sjmisc)

#Solo hematopoyeticas
betas_resultEPIC <- read.csv("betas_result_CellLinesSangre.csv", skip=7)

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

CosmicSSEPIC <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name,]

joinGenes <- read.csv("joinGenes.csv")
joinGenes <- joinGenes$x

CosmicSSEPICGenes <- CosmicSSEPIC[CosmicSSEPIC$GENE_NAME %in% joinGenes,]

sondasYnombres <- DMPs0.15Signi[,c("Name","GencodeBasicV12_NAME","UCSC_RefGene_Name")] #Añadir columan posicion y cromosoma
sondasYnombres$UCSC_RefGene_Name <- sapply(strsplit(sondasYnombres$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
sondasYnombres$GencodeBasicV12_NAME <- sapply(strsplit(sondasYnombres$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])


#CosmicSSEPICNombreGenes <- CosmicSSEPIC$GENE_NAME
#CosmicSSEPICNombreGenes <- unique(CosmicSSEPICNombreGenes)

#CosmicSSEPICSamples <- CosmicCompleteCL$SAMPLE_NAME[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name]
#CosmicSSEPICSamples <- unique(CosmicSSEPICSamples)
#saveRDS(CosmicSSEPICSamples, "CosmicSSEPICSamples.rds")

#for de EPIC.
TCF7L2:TCF4 TRAPPC3L:BET3L
i <- "NLRP3"

pdf("graficas expresion-metilacion/GraficasSangreUCSCUniqInLab.pdf")

#Para ejecutar el for necesitamos cambiar los objetos: listaUCSC, CosmicSS450KGenes y sondasYnombres.
#Leemos una lista de  genes
for(i in joinGenes){
  #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
  if(i %in% CosmicSSEPIC$GENE_NAME){
    
    j <- assign(paste0("gen", i), CosmicSSEPICGenes[CosmicSSEPICGenes$GENE_NAME %in% i,])
    
    sonda <- sondasYnombres[sondasYnombres$UCSC_RefGene_Name %in% i, "Name" ]
    #sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
    
    sonda <- sonda[1]
    
    if(!is_empty(sonda) && sonda  %in% colnames(betas_resultEPIC)){
      expresionGen <- betas_resultEPIC[betas_resultEPIC$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      expresionGen <- as.data.table(expresionGen)
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]
      
      
      
      # plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
      #   geom_point() +
      #   geom_smooth(method='lm') +
      #   stat_cor(method = "spearman" )+ geom_text_repel(label=k$SAMPLE_NAME, max.overlaps = Inf)  + ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" )+ geom_text_repel(label=ifelse(k$SAMPLE_NAME  %in% lineasHematop$Cell.line.name,k$SAMPLE_NAME,""), max.overlaps = Inf)  + ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      
      print(plot)
    }
    #  lista<- append(lista,plot)
    #lista[count] <- plot                               # Store plots in list
    #count <- count + 1
    
  }
}
dev.off()


#lineas Hematopoyeticas:
lineasHematop
lineasHematop <- read.csv("lineasHematop.csv")

lineasHematop$Cell.line.name

expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]

l <- unique(k$SAMPLE_NAME)

outerJoinGenes <- joinGenes[!(joinGenes %in% CosmicSSEPIC$GENE_NAME)]

expresionGen <- as.data.table(expresionGen)





#Bvals binomial:
setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

bVals <- readRDS("bVals.rds")
bVals<- as.data.frame(bVals)
colnames(bVals) <- targets$EPIC.Sample.Name
bValsSS <- bVals[rownames(bVals) %in% c("cg06055229","cg20033583","cg19028462","cg13493526"),]

#bValsSS$`204776850059_R01C01` <- ifelse(bValsSS$`204776850059_R01C01` > 0.3, 1, 0)
bValsSS <- ifelse(bValsSS > 0.3, 1, 0)
bVals <- ifelse(bVals > 0.3, 1, 0)

celulasSS <- c()
celulasSS<- as.data.frame(celulasSS)

cg06055229 <- c(CD4[CD4$sondas %in% "cg06055229","beta"],CD42[CD42$sondas %in% "cg06055229","beta"],CD42[CD42$sondas %in% "cg06055229","beta"],
                CD43[CD43$sondas %in% "cg06055229","beta"],monocito1[monocito1$sondas %in% "cg06055229","beta"],monocito2[monocito2$sondas %in% "cg06055229","beta"],
                monocito3[monocito3$sondas %in% "cg06055229","beta"],CD8[CD8$sondas %in% "cg06055229","beta"],CD82[CD82$sondas %in% "cg06055229","beta"],
                CD83[CD83$sondas %in% "cg06055229","beta"],NeuMale[NeuMale$sondas %in% "cg06055229","beta"],NeuFemale[CD83$NeuFemale %in% "cg06055229","beta"],
                NeuFemale2[NeuFemale2$sondas %in% "cg06055229","beta"])



cg20033583 <- c(CD4[CD4$sondas %in% "cg20033583","beta"],CD42[CD42$sondas %in% "cg20033583","beta"],CD42[CD42$sondas %in% "cg20033583","beta"],
                CD43[CD43$sondas %in% "cg20033583","beta"],monocito1[monocito1$sondas %in% "cg20033583","beta"],monocito2[monocito2$sondas %in% "cg20033583","beta"],
                monocito3[monocito3$sondas %in% "cg20033583","beta"],CD8[CD8$sondas %in% "cg20033583","beta"],CD82[CD82$sondas %in% "cg20033583","beta"],
                CD83[CD83$sondas %in% "cg20033583","beta"],NeuMale[NeuMale$sondas %in% "cg20033583","beta"],NeuFemale[CD83$NeuFemale %in% "cg20033583","beta"],
                NeuFemale2[NeuFemale2$sondas %in% "cg20033583","beta"])




cg20033583 <- c(CD4[CD4$sondas %in% "cg20033583","beta"],CD42[CD42$sondas %in% "cg20033583","beta"],CD42[CD42$sondas %in% "cg20033583","beta"],
                CD43[CD43$sondas %in% "cg20033583","beta"],monocito1[monocito1$sondas %in% "cg20033583","beta"],monocito2[monocito2$sondas %in% "cg20033583","beta"],
                monocito3[monocito3$sondas %in% "cg20033583","beta"],CD8[CD8$sondas %in% "cg20033583","beta"],CD82[CD82$sondas %in% "cg20033583","beta"],
                CD83[CD83$sondas %in% "cg20033583","beta"],NeuMale[NeuMale$sondas %in% "cg20033583","beta"],NeuFemale[CD83$NeuFemale %in% "cg20033583","beta"],
                NeuFemale2[NeuFemale2$sondas %in% "cg20033583","beta"])




setwd("/Users/jparra/Documents/ProyectoCovid/data/blueprint/")
dir()
valoresGenes <- read.csv("valoresGenes.csv")
rownames(valoresGenes) <- valoresGenes$sonda
bValsSS<- as.data.frame(bValsSS)
bValsSS$sondas <- rownames(bValsSS)
celulasSS <- merge(bValsSS,valoresGenes, by.x = "sondas", by.y = "sonda") #Quiza merged que ponga 0 en las que no hacen  merge
#colnames(celulasSS)
celulasSS <- celulasSS[,c("L11-C-Lu", "L13-C-Lu", "L14-C-Lu", "L17-C-Lu", "L3-C-Lu" , "L4-C-Lu","L1-C-Lu", "L15-C-Lu" ,"L5-C-Lu"  ,"L7-C-Lu" , "L8-C-Lu" , "L10-C-Lu" ,"L18-C-Lu",
                          "L24-C-Lu" ,"L25-C-Lu" ,"L19-C-Lu" ,"L22-C-Lu" ,"L23-C-Lu", "L14-P-Lu", "L2-P-Lu",  "L4-P-Lu" , "L5-P-Lu" , "L6-P-Lu" , "L7-P-Lu" , "L8-P-Lu" , "L17-P-Lu",
                          "L18-P-Lu" ,"L15-P-Lu", "L25-P-Lu", "L24-P-Lu" ,"L19-P-Lu" ,"L21-P-Lu", "L22-P-Lu","L23-P-Lu",  "CD4" , "CD42", "CD43",  "monocito1" , 
                          "monocito2" , "monocito3" ,"CD8" ,"CD82" ,"CD83","NeuMale","NeuFemale","NeuFemale2" )]

#rownames(celulasSS) <- c("cg06055229","cg20033583","cg19028462","cg13493526")
heatmap.2(as.matrix(celulasSS),scale="none",trace="none", col=greenred,dendrogram="column")



#cuantas de mis 120 sonsas hay en blueprint
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
celulasYsondas <- read.csv("bigwig/celulasEpiProfile/celulasYsondas.csv")
celulasYsondasSS <- celulasYsondas[celulasYsondas$sondas %in% DMPs0.15Signi$Name,]
celulasYsondasSS$sondas
bValsSS <- bVals[rownames(bVals) %in% celulasYsondasSS$sondas,]
bValsSS <- ifelse(bValsSS > 0.3, 1, 0)
bValsSS<- as.data.frame(bValsSS)
bValsSS$sondas <- rownames(bValsSS)
celulasSS <- merge(bValsSS,celulasYsondasSS, by.x = "sondas", by.y = "sondas") #Quiza merged que ponga 0 en las que no hacen  merge
rownames(celulasSS) <- celulasYsondasSS$sondas
heatmap.2(as.matrix(celulasSS),scale="none",trace="none", col=greenred,dendrogram="column")


bValsSS <- bVals[rownames(bVals) %in% celulasYsondas$sondas,]
celulasSS <- merge(bValsSS,celulasYsondas, by.x = "sondas", by.y = "sondas") #Quiza merged que ponga 0 en las que no hacen  merge

#PCA por Sample_Group
par(mfrow=c(1,2))
plotMDS(celulasSS, gene.selection="common", 
        col=targets$Sample_Group)
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)
write.table(celulasSS, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/celulasSS.csv", sep=",", col.names=TRUE)
celulasSS <- read.csv("bigwig/celulasEpiProfile/celulasSS.csv")



hc = hclust(dist(celulasSS))


#Nombres de genes con biomart:
install.packages("biomartr", dependencies = TRUE)
library(biomartr)
??biomartr
#La linea siguietne funciona al ponerle los dos puntos misteriosamente.
biomart <- biomaRt::useMart("ensembl", host = "https://oct2014.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
filters <-biomaRt::listFilters(biomart)
attributes <- biomaRt::listAttributes(biomart)

m1.biomart <- biomaRt::getBM(filters = "", attributes = c("ensembl_gene_id","entrezgene", "external_gene_name", "hgnc_symbol"), values = "",  mart = biomart)



outerJoinGenes <- joinGenes[!(joinGenes %in% CosmicSSEPIC$GENE_NAME)]
m1.biomart$external_gene_name %in% outerJoinGenes
outer <- m1.biomart[m1.biomart$external_gene_name %in% outerJoinGenes,]
unique(outer$external_gene_name)

CosmicSSEPIC[(CosmicSSEPIC$GENE_NAME %in% "TRAPPC3L"),]

CosmicSSEPIC[(CosmicSSEPIC$GENE_NAME %in% "TCF7L2"),]

CosmicSSEPIC[(CosmicSSEPIC$GENE_NAME %in% "NREP"),]

attributes <- biomaRt::listAttributes(biomart)


m1.biomart[(m1.biomart$ensembl_gene_id %in% "ENSG00000234117"),]
m1.biomart[(m1.biomart$external_gene_name %in% "BET3L"),]
m1.biomart[(m1.biomart$hgnc_symbol %in% "BET3L"),]


getwd()
DiferentesNombresGenes <- read.table(file = "graficas expresion-metilacion/DiferentesNombresGenes.txt", header = TRUE, sep='\t',fill = TRUE)
DiferentesNombresGenes <- read_tsv("/graficas expresion-metilacion/DiferentesNombresGenes.txt")

DiferentesNombresGenes[(DiferentesNombresGenes$Previous.symbol %in% outerJoinGenes),]
DiferentesNombresGenes[(DiferentesNombresGenes$Previous.symbol %in% "FRMD4A"),]
DiferentesNombresGenes[(DiferentesNombresGenes$Alias.symbol %in% "FRMD4A1"),]




sondasYnombres[(sondasYnombres$UCSC_RefGene_Name %in% "TCF4"),]
sondasYnombres[(sondasYnombres$GencodeBasicV12_NAME %in% outerJoinGenes),]
sondasYnombres[(sondasYnombres$UCSC_RefGene_Name %in% outerJoinGenes),]



#MethylResolver:

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")


install.packages("remotes")
remotes::install_github("darneson/MethylResolver")
library(MethylResolver)
bVals <- readRDS("bVals.rds")
rownames(MethylSig)   #MethylSig tiene datos de WB
MethylSigDF <- as.data.frame(MethylSig)
#bValsSS <- bVals[rownames(bVals) %in% rownames(MethylSig),]
bVals <- as.data.frame(bVals)

MethylResolver(methylMix = bVals, methylSig = MethylSig,
               betaPrime = TRUE, outputPath = "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/", outputName = 
                 "MethylResolver", doPar = FALSE, numCores = 1, 
               alpha = seq(0.5,0.9,by = 0.05), absolute = TRUE, purityModel = RFmodel)

MethylResolver(methylMix = bVals, methylSig = MethylSigDF,
               betaPrime = TRUE, outputPath = "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/", outputName = 
                 "MethylResolver4", doPar = FALSE, numCores = 1, 
               alpha = seq(0.5,0.9,by = 0.05), absolute = TRUE, purityModel = RFmodel)

#Cambio absolute
MethylResolver(methylMix = bVals, methylSig = MethylSig,
               betaPrime = TRUE, outputPath = "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/", outputName = 
                 "MethylResolver2", doPar = FALSE, numCores = 1, 
               alpha = seq(0.5,0.9,by = 0.05), absolute = FALSE, purityModel = RFmodel)

class(MethylSig)
#Ejecutamos metilresolver con la matriz de referencia de cibersort.


#library(MethylCIBERSORT)
data("StromalMatrix_V2")


RefPheno <- Stromal_v2.pheno

Stromal_v2 <- as.data.frame(Stromal_v2)

colnames(Stromal_v2) <- Stromal_v2.pheno

Stromal_v2SS <- Stromal_v2[1:500,0:6]

MethylResolver(methylMix = bVals, methylSig = Stromal_v2,
               betaPrime = TRUE, outputPath = "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/", outputName = 
                 "MethylResolver3", doPar = TRUE, numCores = 6, 
               alpha = seq(0.5,0.9,by = 0.05), absolute = TRUE, purityModel = RFmodel)



MethylResolver(methylMix = bVals, outputPath = "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/", outputName = 
                 "MethylResolver5", doPar = FALSE, numCores = 1, 
               alpha = seq(0.5,0.9,by = 0.05), absolute = TRUE, purityModel = RFmodel)





#Depmap


setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")
Depmap <- read.csv("CCLE_expression.csv", TRUE, sep = ",")
joinGenes <- read.csv("joinGenes.csv")
joinGenes <- joinGenes$x


#grepl("REFUND", string, fixed=TRUE)
#DepmapSS2 <- Depmap[,grepl(colnames(Depmap), joinGenes, fixed=TRUE)]
#DepmapSS <- Depmap[,grepl(joinGenes , colnames(Depmap))]

#DepmapSS <- Depmap[grepl(colnames(Depmap), joinGenes, fixed=TRUE),]
#  str_contains(colnames(Depmap),joinGenes)
#DepmapSS3 <- Depmap[,str_contains(joinGenes ,colnames(Depmap))]
#DepmapSS <- Depmap[,str_contains(colnames(Depmap),joinGenes)]





#pdf("graficas expresion-metilacion/GraficasDepmap.pdf")


#for(i in joinGenes){
#Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
if(i %in% colnames(DepmapSS)){
  
  j <- assign(paste0("gen", i), CosmicSSEPICGenes[CosmicSSEPICGenes$GENE_NAME %in% i,])
  
  sonda <- sondasYnombres[sondasYnombres$UCSC_RefGene_Name %in% i, "Name" ]
  #sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
  
  sonda <- sonda[1]
  
  if(!is_empty(sonda) && sonda  %in% colnames(betas_resultEPIC)){
    expresionGen <- betas_resultEPIC[betas_resultEPIC$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
    expresionGen <- as.data.table(expresionGen)
    expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
    
    k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
    
    variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
    #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]
    
    
    
    # plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
    #   geom_point() +
    #   geom_smooth(method='lm') +
    #   stat_cor(method = "spearman" )+ geom_text_repel(label=k$SAMPLE_NAME, max.overlaps = Inf)  + ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
    
    plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
      geom_point() +
      geom_smooth(method='lm') +
      stat_cor(method = "spearman" )+ geom_text_repel(label=ifelse(k$SAMPLE_NAME  %in% lineasHematop$Cell.line.name,k$SAMPLE_NAME,""), max.overlaps = Inf)  + ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
    
    
    print(plot)
  }
  #  lista<- append(lista,plot)
  #lista[count] <- plot                               # Store plots in list
  #count <- count + 1
  
}
#}
dev.off()




#Depmap$Y <- cellLine$displayName[ cellLine$depmapId %in% Depmap$X]

#DepmapSS$Y <- cellLine$displayName[ cellLine$depmapId %in% DepmapSS$X]





#cellLinesNames <- Depmap[,c("X","Y")]


#i <- "IRF2BPL"
#j <- c("MCF2L","X")

betas_result450KSS <- read.csv("betas_result_2021-08-06-13-40.csv", skip=7) #450K de pulmon

betas_result450KSS$Sample.Name <- str_replace(betas_result450KSS$Sample.Name, "-", "")


betas_resultEPIC <- read.csv("betas_result_CellLinesSangre.csv", skip=7)

betas_resultEPIC$Sample.Name <- str_replace(betas_resultEPIC$Sample.Name, "-", "")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
sondasYnombres <- DMPs0.15Signi[,c("Name","GencodeBasicV12_NAME","UCSC_RefGene_Name")] #Añadir columan posicion y cromosoma
sondasYnombres$UCSC_RefGene_Name <- sapply(strsplit(sondasYnombres$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
sondasYnombres$GencodeBasicV12_NAME <- sapply(strsplit(sondasYnombres$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])

library(stringr)
library(limma)   #Zscore
library(sjmisc)  # is_empty
library(ggpubr)  #ggplot #stat_cor
library(ggrepel)  #geom_text_repel
library(data.table)

#rm(Depmap5)
#rm(sonda)



colnames(Depmap)<- sapply(strsplit(colnames(Depmap), '.', fixed=T), function(x) x[1])
DepmapSS<-   Depmap[,colnames(Depmap) %in% joinGenes]
#DepmapSS$X<-   Depmap$X[colnames(Depmap) %in% joinGenes]
DepmapSS$X<-   Depmap$X
#colnames(DepmapSS)<- colnames2 

cellLine <- read.csv("cell-line-selector.csv", TRUE, sep = ",")
DepmapIdToCellLine <- cellLine[ cellLine$depmapId %in% Depmap$X,c("displayName","depmapId")]


pdf("graficas expresion-metilacion/GraficasDepmapGencodeAllLung.pdf")
#i <- "EBF1"
for(i in joinGenes){
  if(i %in% colnames(DepmapSS)){
    
    #Depmap5 <- DepmapSS[,colnames(DepmapSS) %in% c(i,"X")]
    
    Depmap5 <- DepmapSS[,colnames(DepmapSS) %in% "X"]
    Depmap5 <- as.data.frame(Depmap5)
    
    Depmap5$gen <- DepmapSS[,colnames(DepmapSS) %in% i]
    
    #Debo pasar los datos de expresion de TPM(Transcritos por millon) a Z_SCORE(Normalizados): library(limma)
    
    Depmap5$gen <- zscore(Depmap5$gen, dist="gamma", shape=0.5)
    
    
    k <- merge(DepmapIdToCellLine, Depmap5, by.x = "depmapId", by.y = "Depmap5")
    
    #sonda <- sondasYnombres[sondasYnombres$UCSC_RefGene_Name %in% i, "Name" ]
    sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
    
    sonda <- sonda[1]
    
    if(!is_empty(sonda) && sonda  %in% colnames(betas_result450KSS)){
      
      #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeCompV12_NAME),"GencodeCompV12_Group"]
      
      expresionGen <- betas_result450KSS[betas_result450KSS$Sample.Name %in% k$displayName ,c("Sample.Name", sonda)]
      expresionGen <- as.data.table(expresionGen)
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      l <- merge(k, expresionGen, by.x = "displayName", by.y = "Sample.Name")
      
      plot <- ggplot(l,aes(x = l[ ,sonda],y = l[ ,"gen"])) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        #geom_text_repel(label=l$displayName, max.overlaps = Inf) +
        ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      print(plot)
    }
  }
}
dev.off()


plot <- ggplot(l,aes(x = l[ ,sonda],y = l[ ,"gen"])) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman" ) +
  geom_text_repel(label=ifelse(l$displayName  %in% lineasHematop$Cell.line.name,l$displayName,""), max.overlaps = Inf) +
  ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))

print(plot)


dim(k[k$displayName %in% expresionGen$Sample.Name,])
dim(betas_resultEPIC$Sample.Name %in% k$displayName )





#Ver la diferencia de metilacion del gen que es marcador de los fibroblastos patologicos
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

ann850kgen <- ann850k$Name[ann850k$UCSC_RefGene_Name %in%  "CTHRC1"]
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
DMPs0.15SigniGen <- DMPs0.15Signi[DMPs0.15Signi$Name %in% ann850kgen]

bVals <- readRDS("bVals.rds")
bValsDF <- as.data.frame(bVals)


BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")

BetaDiferenciaMediaSS <- BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% rownames(DMPs0.15Signi),]

BetaDiferenciaMediaGen <- BetaDiferenciaMedia$Diferencia[rownames(BetaDiferenciaMedia) %in% ann850kgen]


BetaDiferenciaMediaGen <- BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% rownames(DMPs0.15Signi),]
write.table(BetaDiferenciaMediaGen, file="BetaDiferenciaMediaCpGs.csv", sep=",", row.names=TRUE)


bValsDFGen <- bValsDF[rownames(bValsDF) %in% ann850kgen,]

#heatmap.2 gplots
heatmap.2(as.matrix(bValsDFGen),scale="none",trace="none", col=greenred,dendrogram="column")





#Graficas de metilacion expresion para todas las lineas celulares de EPIC en methylationDB.
library(tidyverse) #read_tsv
library(data.table) #as.data.table
library(ggpubr)  #ggplot #stat_cor

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

betas_resultEPIC <- read.csv("betas_resultAllCellLines.csv", skip=7)


joinGenes <- read.csv("joinGenes.csv")
joinGenes <- joinGenes$x


CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name,]

CosmicSSEPICGenes <- CosmicSSAll[CosmicSSAll$GENE_NAME %in% joinGenes,]

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
sondasYnombres <- DMPs0.15Signi[,c("Name","GencodeBasicV12_NAME","UCSC_RefGene_Name")] #Añadir columan posicion y cromosoma
sondasYnombres$UCSC_RefGene_Name <- sapply(strsplit(sondasYnombres$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
sondasYnombres$GencodeBasicV12_NAME <- sapply(strsplit(sondasYnombres$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])

DMPs0.15Signi$UCSC_RefGene_Name <- sapply(strsplit(sondasYnombres$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_NAME <- sapply(strsplit(sondasYnombres$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])

DMPs0.15Signi$UCSC_RefGene_Group <- sapply(strsplit(sondasYnombres$UCSC_RefGene_Group, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_Group <- sapply(strsplit(sondasYnombres$GencodeBasicV12_Group, ';', fixed=T), function(x) x[1])

i <- "A2ML1"

pdf("graficas expresion-metilacion/GraficaAllCellLinesUCSCOur.pdf")

#Para ejecutar el for necesitamos cambiar los objetos: listaUCSC, CosmicSS450KGenes y sondasYnombres.
#Leemos una lista de  genes
for(i in joinGenes){
  #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
  if(i %in% CosmicSSAll$GENE_NAME){
    
    j <- assign(paste0("gen", i), CosmicSSEPICGenes[CosmicSSEPICGenes$GENE_NAME %in% i,])
    
    sonda <- sondasYnombres[sondasYnombres$UCSC_RefGene_Name %in% i, "Name" ]
    #sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
    
    sonda <- sonda[1]
    
    if(!is_empty(sonda) && sonda  %in% colnames(betas_resultEPIC)){
      expresionGen <- betas_resultEPIC[betas_resultEPIC$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      expresionGen <- as.data.table(expresionGen)
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]
      
      
      
      # plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
      #   geom_point() +
      #   geom_smooth(method='lm') +
      #   stat_cor(method = "spearman" )+ geom_text_repel(label=k$SAMPLE_NAME, max.overlaps = Inf)  + ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
        labs(x=sonda) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        # geom_text_repel(label=k$SAMPLE_NAME) +
        # geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab$Cell.line.name,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
        geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% lineasHematop$Cell.line.name,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
        ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      
      print(plot)
    }
    #  lista<- append(lista,plot)
    #lista[count] <- plot                               # Store plots in list
    #count <- count + 1
    
  }
}
dev.off()

geom_text(label=k$SAMPLE_NAME)
geom_text(aes(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab,k$SAMPLE_NAME,"")),hjust = 1.1)

CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% joinGenes,]




betas_result450KSS <- read.csv("betas_result_2021-08-06-13-40.csv", skip=7) #450K de pulmon

betas_result450KSS$Sample.Name <- str_replace(betas_result450KSS$Sample.Name, "-", "")


CosmicSS450K <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_result450KSS$Sample.Name,]
CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name,]



pdf("graficas expresion-metilacion/GraficaLung450GencodeInLab.pdf")

#Para ejecutar el for necesitamos cambiar los objetos: listaUCSC, CosmicSS450KGenes y sondasYnombres.
#Leemos una lista de  genes
for(i in joinGenes){
  #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
  if(i %in% CosmicSS450KGenes$GENE_NAME){
    
    j <- assign(paste0("gen", i), CosmicSS450KGenes[CosmicSS450KGenes$GENE_NAME %in% i,])
    
    #sonda <- sondasYnombres[sondasYnombres$UCSC_RefGene_Name %in% i, "Name" ]
    sonda <- sondasYnombres[sondasYnombres$GencodeBasicV12_NAME %in% i, "Name" ]
    
    sonda <- sonda[1]
    
    if(!is_empty(sonda) && sonda  %in% colnames(betas_result450KSS)){
      expresionGen <- betas_result450KSS[betas_result450KSS$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      expresionGen <- as.data.table(expresionGen)
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]
      
      
      plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
        labs(x=sonda) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        # geom_text_repel(label=k$SAMPLE_NAME) +
        geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
        #geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% lineasHematop$Cell.line.name,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
        ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))
      
      
      print(plot)
    }
    
  }
}
dev.off()


#Uniendo pulmon(450K) con hematopeyeticas(EPIC)
setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

library(tidyverse) #read_tsv
library(data.table) #as.data.table
library(ggpubr)  #ggplot #stat_cor
library(ggrepel)  #geom_text_repel

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

betas_resultEPIC <- read.csv("betas_result_CellLinesSangre.csv", skip=7) #Epic de sangre
betas_result450KSS <- read.csv("betas_result_2021-08-06-13-40.csv", skip=7) #450K de pulmon

betasSS1 <- betas_resultEPIC[,c("Sample.Name","cg13493526")]
betasSS2 <- betas_result450KSS[,c("Sample.Name","cg13493526")]

betasAll <- rbind(betasSS1,betasSS2)
CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]

CosmicSSAll <- CosmicSSAll[CosmicSSAll$GENE_NAME%in% "IRF2BPL", ]

betasAll <- as.data.table(betasAll)
betasAll <- betasAll[,lapply(.SD, mean), by = Sample.Name]

k <- merge(CosmicSSAll,betasAll , by.x = "SAMPLE_NAME", by.y = "Sample.Name")

plot <- ggplot(k,aes(x = cg13493526,y = Z_SCORE)) +
  labs(x="cg13493526") +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman" ) +
  geom_text_repel(label=k$SAMPLE_NAME) +
  #geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
  #geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% lineasHematop$Cell.line.name,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
  ggtitle(paste("Grafica expresion-metilacion para el gen:IRF2BPL"))


print(plot)







#Uniendo todo 450K con todo EPIC
library(tidyverse) #read_tsv

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

betas_result450K<- read.csv("betas_result_450KAll1433.csv", skip=7) #Todo 450K

betas_resultEPIC <- read.csv("betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC

#Cogemos los datos de epic para las sonas de 450K
betas_resultEPICSS <- betas_resultEPIC[colnames(betas_resultEPIC) %in% colnames(betas_result450K)]

betasAll <- rbind(betas_result450K,betas_resultEPICSS)

i <- "NLRP3"

#NLRP3 cg18126557
betasAll <- betasAll[,c("Sample.Name","cg18126557")]

CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]

CosmicSSAll <- CosmicSSAll[CosmicSSAll$GENE_NAME %in% "NLRP3", ]

library(data.table) #as.data.table
betasAll <- as.data.table(betasAll)
betasAll <- betasAll[,lapply(.SD, mean), by = Sample.Name]

k <- merge(CosmicSSAll,betasAll , by.x = "SAMPLE_NAME", by.y = "Sample.Name")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
DMPs0.15Signi$UCSC_RefGene_Name <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_NAME <- sapply(strsplit(sondasYnombres$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])

DMPs0.15Signi$UCSC_RefGene_Group <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Group, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_Group <- sapply(strsplit(sondasYnombres$GencodeBasicV12_Group, ';', fixed=T), function(x) x[1])

variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
#variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]

library(data.table) #as.data.table
library(ggpubr)  #ggplot #stat_cor
library(ggrepel)  #geom_text_repel
plot <- ggplot(k,aes(x = cg18126557,y = Z_SCORE)) +
  labs(x="cg18126557") +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor(method = "spearman" ) +
  geom_text_repel(label=k$SAMPLE_NAME) +
  
  #geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% lineasHematop$Cell.line.name,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
  #ggtitle(paste("Grafica expresion-metilacion para el gen:NLRP3"))
  ggtitle(paste("Grafica expresion-metilacion para el gen:",i,"que se encuenta en :",variablePromotor))

print(plot)



#Para todos los genes

library(tidyverse) #read_tsv

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

betas_result450K<- read.csv("betas_result_450KAll1433.csv", skip=7) #Todo 450K

betas_resultEPIC <- read.csv("betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
betas_resultEPICSS <- betas_resultEPIC[colnames(betas_resultEPIC) %in% colnames(betas_result450K)]

betasAll <- rbind(betas_result450K,betas_resultEPICSS)
CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]
#k <- merge(CosmicSSAll,betasAll , by.x = "SAMPLE_NAME", by.y = "Sample.Name")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
DMPs0.15Signi$UCSC_RefGene_Name <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_NAME <- sapply(strsplit(DMPs0.15Signi$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])

DMPs0.15Signi$UCSC_RefGene_Group <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Group, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_Group <- sapply(strsplit(DMPs0.15Signi$GencodeBasicV12_Group, ';', fixed=T), function(x) x[1])



joinGenes <- read.csv("joinGenes.csv")
joinGenes <- joinGenes$x

#CosmicSS450K <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_result450K$Sample.Name,]
#CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% joinGenes,]

#betasAll[,colnames(betasAll) %in% "cg26818401"]
#betas_resultEPIC[,colnames(betas_resultEPIC) %in% "cg26818401"]
#betas_result450K[,colnames(betas_result450K) %in% "cg26818401"]


pdf("graficas expresion-metilacion/GraficaAllUCSC.pdf")
library(ggpubr)  #ggplot #stat_cor
library(ggrepel)  #geom_text_repel
i <- "ZNF608" 
for(i in joinGenes){
  #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
  if(i %in% CosmicSSAll$GENE_NAME){
    
    
    
    #j <- assign(paste0("gen", i), CosmicSS450KGenes[CosmicSS450KGenes$GENE_NAME %in% i,])
    j <- assign(paste0("gen", i), CosmicSSAll[CosmicSSAll$GENE_NAME %in% i,])
    
    sonda <- DMPs0.15Signi[DMPs0.15Signi$UCSC_RefGene_Name %in% i, "Name" ]
    #sonda <- DMPs0.15Signi[DMPs0.15Signi$GencodeBasicV12_NAME %in% i, "Name" ]
    
    sonda <- sonda[1]
    
    if(!is_empty(sonda) && sonda  %in% colnames(betasAll)){
      expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      expresionGen <- as.data.table(expresionGen)
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]
      
      plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
        labs(x=sonda) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        geom_text_repel(label=k$SAMPLE_NAME) +
        #geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
        ggtitle(paste("Expression-methylation graph for the gene:",i,"that is found in:",variablePromotor))
      
      print(plot)
      
      
    }
  } #a
}  

dev.off()








#Pruebas scale

x1 <- c(0.01,0.02,0.03)
sd1 <- sd(x1)

x2 <- c(0.80,0.97,0.99)
sd2 <- sd(x2)


#sondas en 450K
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sondas850K <- rownames(ann850k)
sondas <- c("cg07465999","cg24246165","cg05648629","cg11787526","cg24679427","cg10404427","cg21124310","cg24705426","cg17296589","cg07525294","cg16394018","cg09922561","cg22483613","cg15367056","cg22980156","cg20775850","cg06725421","cg08423585","cg08459290","cg27127608","cg09606770","cg24480152","cg03608224","cg21741010","cg14547726","cg10743544","cg17920246","cg00524377","cg14253027","cg24874003","cg24997442","cg10996039","cg16325432","cg06360465","cg02101707","cg13493526","cg25452407","cg21413173","cg03398002","cg12251779","cg08640619","cg22715761","cg07416291","cg19028462","cg25113767","cg18126557","cg11895835","cg23256675","cg15745560","cg13980607","cg05122159","cg24598922","cg05291674","cg09485632","cg26818401","cg10985150","cg05490233","cg27365602","cg24822602","cg25855753","cg20033583","cg02717112","cg02531640","cg05472536","cg08884752","cg27057797","cg11189357","cg02471028","cg20017590","cg20263678","cg14568422","cg26554902","cg00814721","cg24388610","cg03719830","cg19870403","cg19902005","cg27496650","cg10298187","cg16811988","cg04309595","cg00285941","cg10608681","cg09025210","cg25458175","cg18521743","cg06055229")
sondas850K$sondas850K[sondas850K$sondas850K %in% sondas]

sondas[sondas %in% sondas450K]
sondas <- as.data.frame(sondas)
sondas850K <- as.data.frame(sondas850K)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sondas450K <- rownames(ann450k)


genesCoding <- c("A2ML1","ACYP2","ADCY9","ADORA3","ANKRD55","ATG7","ATP11A","BET3L","BIRC6","CALCB","CAV1","CD38","CDK6",
                 "CDRT4","CHST11","CTNND1","DHRS12","DSTN","7SK","DTNB","EBF1","EIF4E3","EMP2","FAM101A",
                 "FAM26E","FMNL3","FRMD4A","FSTL1","GFRA2","GNG7","GPX3","HOXA3","HOXA3","HYAL1","IPO9",
                 "IRF2BPL","ITPR1","KAZN","KIAA0930","KIRREL","LDB3","MAST4","MCF2L","MED12L","NCOR2",
                 "NLRP3","NREP","NUPR1","P2RY14","PDE1C","PHACTR2","PIAS1","PIK3R1","PIP5K1B","PLEKHG1",
                 "PLXND1","PPARGC1A","PREP","PRR34","RAB11FIP5","RAPGEF2","RBPMS","RN5S184","SAP30BP",
                 "SKI","SLC24A3","SLC45A1","SPON2","SPTBN1","TBC1D10A","TCF4","TEX30","THADA","TM9SF4",
                 "TMEM132B","TMIGD3","TMX2-CTNND1","TOX","TPPP3","TRAK1","TRAM2","UNC80","WDR70","ZBTB38",
                 "ZBTB45","ZNF608")

genesMEMORY_CD4_TCELL <- c("ACP5","ADD1","ADK","ADO","AGPAT5","ALDH2","ALOX5","ARHGAP17","ARHGEF7","BACH2",
                           "BANK1","BASP1","BCL11A","BHLHE41","BLNK","BMP2K","BTN2A2","C5orf15","CA2",
                           "CACNA1A","CAT","CCNB1IP1","CCT2","CD180","CD1D","CD200","CD38","CD72","CD74",
                           "CD79B","CDK14","CEPT1","CHFR","CKAP4","CLEC4A","COBLL1","COL9A3","COQ2","COX15",
                           "COX17","CSNK2A1","CTNNA1","CTSH","CYBB","DGKD","DLAT","DUS2","E2F5","EAF2",
                           "EIF2AK3","EIF2B1","EMC8","ERCC1","FAM30A","FCER2","FCGR2C","GABARAP","GM2A",
                           "GNA12","GNG7","GSN","GUSBP11","H2AC6","H2BC12","H2BC9","HESX1","HHEX","HLA-DMA",
                           "HLA-DMB","HLA-DOB","HLA-DQB1","HPS5","HS3ST1","HSD17B4","IGHA1","IGHD","IGHM",
                           "IGKC","IGKV1D-13","IGKV4-1","IGLL3P","IL13RA1","IL4R","ING1","IRAK1","IRF4",
                           "IRF8","ISG20","ITPR1","JCHAIN","JUN","KCNN4","KIAA0040","LAMP5","LAT2","LMO2",
                           "LY86","LYN","LYST","MAP4K2","MAPK12","MCM5","MED14","MEF2C","METTL7A","MICAL3",
                           "MICALL1","MRPL15","MRPL40","MS4A1","MZB1","NADK","NASP","NT5E","NUMB","OAS1",
                           "OPN3","ORC5","OSBPL10","P2RX1","P2RX5","PAPSS1","PARM1","PARP1","PAX5","PDLIM1",
                           "PECAM1","PIP5K1B","PKIG","PLCG2","POU2AF1","PSEN2","PSMB6","PTK2B","PTPN6","PUS1",
                           "QRSL1","RAB31","RABGEF1","RASGRP3","REXO4","RHOB","RUBCN","RUFY1","S100A8","SAV1",
                           "SCN3A","SCO2","SCPEP1","SEC23B","SEL1L3","SERPINF1","SHMT2","SIDT1","SIDT2",
                           "SIPA1L3","SLC15A2","SLC15A3","SLC37A1","SNRNP25","SP110","SPIB","SQLE","STAP1",
                           "STX7","SWAP70","SYK","SYPL1","TBC1D9","TCF4","TERF2","TFEB","TLE1","TLR7","TMEM14B",
                           "TMEM156","TMEM243","TMEM62","TMEM70","TMUB2","TNFRSF17","TOR3A","TRAK1","TSPAN13",
                           "TSPAN3","UBE2J1","UROD","USP22","WARS1","XIST","YBX3","ZBTB5","ZDHHC14","ZNF106",
                           "ZNF165","ZNF232","ZNF318")

genesTIMULATED_IRF4_KO_BCELL <- c("ABLIM1","ALDH2","ALPK2","AMZ2P1","ANKRD11","ANKRD50","ARHGEF2","ASMTL",
                                  "ATF7IP","AURKA","BANP","BRD1","C1QL1","CACNA2D2","CAV1","CD40","CD70",
                                  "CD86","CDC42EP4","CDC42SE1","CDK2","CDYL","CENATAC","CEP350","CEP68",
                                  "CETP","CGN","CHD2","CLCF1","CLDND1","CLIP2","CREB5","CREM","CXCR4",
                                  "CYTH1","DCTN2","DMXL1","DNAJC1","DPP9","DSG2","DSTN","DSTNP2","DUSP4",
                                  "DYNLRB2","EBF1","EMP2","ENSG00000268472","FAM107B","FBH1","FBN2",
                                  "FBXO7","FLJ13224","FSTL3","FXR1","GCNA","GLT8D1","GNA12","GNG4",
                                  "GSTA4","GTF2IRD2","HES1","HMGB3","HS1BP3","HTR2B","HUS1","IFRD1",
                                  "IGSF8","IL12B","IL18R1","IL23A","ITPR1","JOSD1","KCNE5","KIF13A",
                                  "KIF21A","KIF2C","KLF13","KPNA1","KRBOX4","KSR1","LINC00294",
                                  "LINC00410","LOC440896","LRIG1","LUZP1","LYSMD2","MATN2","MDFIC",
                                  "MED13","MED29","MIR503HG","MLLT11","MMP26","MTRFR","MVB12B","NCOR2",
                                  "NEFH","NKAPD1","NR4A1","NXPH3","OR3A3","OTUD7B","PABPC3","PAIP2",
                                  "PBX2","PCDHB17P","PDE4B","PDE4D","PDK4","PDRG1","PHC3","PIK3R3",
                                  "PIP5K1B","PITPNC1","PKIA","PLXNC1","POU2F1","PPM1A","PRKCE","PRRC2B",
                                  "PTPN1","RAB3IP","RAB5B","RAPGEF6","RARRES1","RASA1","RBBP6","RBM23",
                                  "RETREG1","RFC2","RGS2","RHOB","RNF168","S100A3","SAP130","SEC22A",
                                  "SEMA3F","SEMA7A","SH3PXD2A","SKAP1","SLC12A4","SLC5A6","SNORA71B",
                                  "SORL1","SPINDOC","SRSF12","SSH1","ST20-AS1","ST3GAL1","ST6GALNAC6",
                                  "ST8SIA6-AS1","SV2B","SYNPO","TAF6","TAMALIN","TBL1X","TFAP2D","TGIF1",
                                  "THAP3","TKTL1","TLCD3A","TMEM145","TMOD2","TMPRSS4","TNFAIP1","TNPO2",
                                  "TP53BP1","TPRA1","TPST1","TRABD2A","TRIM15","TSHZ1","TSPAN13","TSPAN2",
                                  "TTK","TTYH2","UBAC2","UNKL","USP12","USP54","WHRN","WNK1","WWP1","XCL1",
                                  "ZADH2","ZC2HC1A","ZC3H4","ZMIZ1","ZNF250","ZNF449","ZNF516","ZNF598",
                                  "ZNF627","ZNF79","ZSCAN12")

unique(genesCoding)

genesMEMORY_CD4_TCELL[genesMEMORY_CD4_TCELL %in% genesCoding]

genesTIMULATED_IRF4_KO_BCELL[genesTIMULATED_IRF4_KO_BCELL %in% genesCoding]


#wilcoxon:


#comprobar normalidad.
shapiro.test(targets$AGE) # p-value = 0.1694. Si es normal.

#Utilizo Wilcoxon para comparar dos variables independientes, una categórica y otra numerica(con distribucion normal) con dos categorías.
wilcox.test(targets$AGE ~ targets$TISSUE) 
#p-value = 0.9862. No hay relacion

t.test(targets$AGE ~ targets$TISSUE, data = targets) #p-value = 0.9597. Con t.test tampoco hay relacion, prefiero wilcoxon porque t.test es para variables dependientes(por ejemplo antes y despues)¿?


kruskal.test(targets$GENDER ~ targets$TISSUE, data = targets)
#p-value = 0.5574. No hay relacion


#Ward2
data <- c(14,11,4,5)
data <- matrix(data, ncol=2)
fisher.test(data)


data <- c(10,0,6,18)
data <- matrix(data, ncol=2)
fisher.test(data)



#Control vs pneumonia

dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")
bVals <- readRDS("bVals.rds")
colnames(bVals) <- targets$EPIC.Sample.Name

MuestrasCP <- targets$EPIC.Sample.Name[targets$TISSUE %in% "CONTROL"] 
MuestrasCP <- c(MuestrasCP,"L8-P-Lu","L6-P-Lu","L22-P-Lu")

bValsCP <- bVals[,colnames(bVals) %in% MuestrasCP]

targetsCP <- targets[targets$EPIC.Sample.Name %in% c("L11-C-Lu","L13-C-Lu","L14-C-Lu","L17-C-Lu","L3-C-Lu","L4-C-Lu","L1-C-Lu","L15-C-Lu",
                                                     "L5-C-Lu","L7-C-Lu","L8-C-Lu","L10-C-Lu","L18-C-Lu","L24-C-Lu","L25-C-Lu","L19-C-Lu",
                                                     "L22-C-Lu","L23-C-Lu","L8-P-Lu","L6-P-Lu","L22-P-Lu"),]

# this is the factor of interest
tissue <- factor(targetsCP$TISSUE)
# this is the individual effect that we need to account for 
#gender <- factor(targets$GENDER) 
# use the above to create a design matrix
design <- model.matrix(~0+tissue, data=targetsCP) 
colnames(design) <- levels(tissue)
# fit the linear model 
#dim(mVals)
library(limma)
fit <- lmFit(bValsCP, design)
# create a contrast matrix for specific comparisons                     
contMatrix <- makeContrasts(CONTROL-COVID,
                            levels=design) 
contMatrix




fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))


#Este match es muy imprtante para que anote en orden
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850kSub <- ann850k[match(rownames(bVals),ann850k$Name),]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann850kSub)
#DMPsSigni <- subset(DMPs, adj.P.Val < 0.05)
DMPsSigni <- DMPs[DMPs$adj.P.Val < 0.05,]
dim(DMPsSigni)
#De todas las sondas filtramos las que tengan un p.valor ajustado menor de 0.05,
#osea  las que tengan una correlacion segun la regresion lineal de limma. 

bValsCP <- bVals[,colnames(bVals) %in% MuestrasCP]

colnames(bValsCP)

suma_renglones <- function (x) {
  y=integer(nrow(x))
  z=integer(nrow(x))
  r=integer(nrow(x))
  names <- colnames(x)
  for(i in 1:nrow(x)){
    for (j in 1:ncol(x)) {
      if(is.numeric(x[i,j]) && (names[j]=="L1-C-Lu"
                                || names[j]=="L3-C-Lu"
                                || names[j]=="L4-C-Lu"
                                || names[j]=="L5-C-Lu"
                                || names[j]=="L7-C-Lu"
                                || names[j]=="L8-C-Lu"
                                || names[j]=="L10-C-Lu"
                                || names[j]=="L11-C-Lu"
                                || names[j]=="L13-C-Lu"
                                || names[j]=="L14-C-Lu"
                                || names[j]=="L15-C-Lu"
                                || names[j]=="L17-C-Lu"
                                || names[j]=="L18-C-Lu"
                                || names[j]=="L19-C-Lu"
                                || names[j]=="L22-C-Lu"
                                || names[j]=="L23-C-Lu"
                                || names[j]=="L24-C-Lu"
                                || names[j]=="L25-C-Lu" )) 
      {
        y[i] = y[i]+(x[i,j])
      }else if (is.numeric(x[i,j]) && (names[j]=="L6-P-Lu"
                                       || names[j]=="L8-P-Lu"
                                       || names[j]=="L22-P-Lu")) {
        z[i] = z[i]+(x[i,j])
      }
    }
    dim(y)
    dim(x)
    y[i] <-y[i]/18
    z[i] <-z[i]/3
    r[i] = y[i]-z[i]
  }
  
  x = as.data.frame(x)
  x[, "mediaControl"] <- y
  x[, "mediaCovid"] <- z
  x[, "Diferencia"] <- r
  
  return(x)
}



BetaDiferenciaMedia = suma_renglones(bValsCP)

#BetaDiferenciaMediaSS <- subset(BetaDiferenciaMedia, Diferencia > 0.25 | Diferencia < -0.25)
BetaDiferenciaMediaSS <- BetaDiferenciaMedia[(BetaDiferenciaMedia$Diferencia > 0.25 | BetaDiferenciaMedia$Diferencia < -0.25),]

dim(BetaDiferenciaMediaSS)
#Me da 11200 sondas con diferencia entre medias mayor de 0.15
#saveRDS(BetaDiferenciaMediaSS, "BetaDiferenciaMediaSSCP.rds")
#BetaDiferenciaMediaSS <- readRDS("BetaDiferenciaMediaSSCP.rds")

DMPsSSCP <- DMPsSigni[match(rownames(BetaDiferenciaMediaSS),rownames(DMPsSigni)),]

#DMPsSSCP <- DMPsSigni[match(rownames(DMPsSigni),rownames(BetaDiferenciaMediaSS)),]

DMPsSSCP <-DMPsSSCP[!is.na(DMPsSSCP$chr),]
dim(DMPsSSCP)


DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
##DMPsConVsPneu <- DMPs0.15Signi[rownames(DMPs0.15Signi) %in% rownames(BetaDiferenciaMediaSS),]

UCSC<- strsplit(as.character(DMPsSSCP$UCSC_RefGene_Name), ";")

UCSC <-unlist(UCSC)

Gencode<- strsplit(as.character(DMPsSSCP$GencodeBasicV12_NAME), ";")

Gencode <-unlist(Gencode)

genes <- c(UCSC,Gencode)

genes <- unique(genes)

write.table(genes, file="genes.csv", sep=",", row.names=FALSE)



BetaDiferenciaMediaSSCP <- BetaDiferenciaMediaSS[match(rownames(DMPsSSCP),rownames(BetaDiferenciaMediaSS)),]
write.table(BetaDiferenciaMediaSSCP, file="BetaDiferenciaMediaSSCP.csv", sep=",", row.names=TRUE)




#Control vs acute DAD
#acute DAD: 4,5,14,15,17,18,19,23,24,25

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")
bVals <- readRDS("bVals.rds")
colnames(bVals) <- targets$EPIC.Sample.Name

MuestrasCP <- targets$EPIC.Sample.Name[targets$TISSUE %in% "CONTROL"] 
MuestrasCP <- c(MuestrasCP,"L4-P-Lu","L5-P-Lu","L14-P-Lu","L15-P-Lu","L17-P-Lu","L18-P-Lu","L19-P-Lu","L23-P-Lu","L24-P-Lu","L25-P-Lu")

bValsCP <- bVals[,colnames(bVals) %in% MuestrasCP]

targetsCP <- targets[targets$EPIC.Sample.Name %in% c("L11-C-Lu","L13-C-Lu","L14-C-Lu","L17-C-Lu","L3-C-Lu","L4-C-Lu","L1-C-Lu","L15-C-Lu",
                                                     "L5-C-Lu","L7-C-Lu","L8-C-Lu","L10-C-Lu","L18-C-Lu","L24-C-Lu","L25-C-Lu","L19-C-Lu",
                                                     "L22-C-Lu","L23-C-Lu","L4-P-Lu","L5-P-Lu","L14-P-Lu","L15-P-Lu","L17-P-Lu","L18-P-Lu",
                                                     "L19-P-Lu","L23-P-Lu","L24-P-Lu","L25-P-Lu"),]

# this is the factor of interest
tissue <- factor(targetsCP$TISSUE)
# this is the individual effect that we need to account for 
#gender <- factor(targets$GENDER) 
# use the above to create a design matrix
design <- model.matrix(~0+tissue, data=targetsCP) 
colnames(design) <- levels(tissue)
# fit the linear model 
#dim(mVals)
library(limma)
fit <- lmFit(bValsCP, design)
# create a contrast matrix for specific comparisons                     
contMatrix <- makeContrasts(CONTROL-COVID,
                            levels=design) 
contMatrix




fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))


#Este match es muy imprtante para que anote en orden
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850kSub <- ann850k[match(rownames(bVals),ann850k$Name),]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann850kSub)


#bValsCP <- bVals[,colnames(bVals) %in% MuestrasCP]

colnames(bValsCP)

suma_renglones <- function (x) {
  y=integer(nrow(x))
  z=integer(nrow(x))
  r=integer(nrow(x))
  names <- colnames(x)
  for(i in 1:nrow(x)){
    for (j in 1:ncol(x)) {
      if(is.numeric(x[i,j]) && (names[j]=="L1-C-Lu"
                                || names[j]=="L3-C-Lu"
                                || names[j]=="L4-C-Lu"
                                || names[j]=="L5-C-Lu"
                                || names[j]=="L7-C-Lu"
                                || names[j]=="L8-C-Lu"
                                || names[j]=="L10-C-Lu"
                                || names[j]=="L11-C-Lu"
                                || names[j]=="L13-C-Lu"
                                || names[j]=="L14-C-Lu"
                                || names[j]=="L15-C-Lu"
                                || names[j]=="L17-C-Lu"
                                || names[j]=="L18-C-Lu"
                                || names[j]=="L19-C-Lu"
                                || names[j]=="L22-C-Lu"
                                || names[j]=="L23-C-Lu"
                                || names[j]=="L24-C-Lu"
                                || names[j]=="L25-C-Lu" )) 
      {
        y[i] = y[i]+(x[i,j])
      }else if (is.numeric(x[i,j]) && (names[j]=="L4-P-Lu"
                                       || names[j]=="L5-P-Lu"
                                       || names[j]=="L14-P-Lu"
                                       || names[j]=="L15-P-Lu"
                                       || names[j]=="L17-P-Lu"
                                       || names[j]=="L18-P-Lu"
                                       || names[j]=="L19-P-Lu"
                                       || names[j]=="L23-P-Lu"
                                       || names[j]=="L24-P-Lu"
                                       || names[j]=="L25-P-Lu")) {
        z[i] = z[i]+(x[i,j])
      }
      
    }
    dim(y)
    dim(x)
    y[i] <-y[i]/18
    z[i] <-z[i]/10
    r[i] = y[i]-z[i]
  }
  
  x = as.data.frame(x)
  x[, "mediaControl"] <- y
  x[, "mediaCovid"] <- z
  x[, "Diferencia"] <- r
  
  return(x)
}



BetaDiferenciaMedia = suma_renglones(bValsCP)

BetaDiferenciaMediaSS <- subset(BetaDiferenciaMedia, Diferencia > 0.15 | Diferencia < -0.15)
dim(BetaDiferenciaMediaSS)
#Me da 265 sondas con diferencia entre medias mayor de 0.15
saveRDS(BetaDiferenciaMediaSS, "BetaDiferenciaMediaSSCP.rds")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
DMPsConVsAcuteDAD <- DMPs0.15Signi[rownames(DMPs0.15Signi) %in% rownames(BetaDiferenciaMediaSS),]
#En comun 30





DMPsSSCP <- DMPs[match(rownames(BetaDiferenciaMediaSS),rownames(DMPs)),]
dim(DMPsSSCP)

DMPsSSCPSigni <- subset(DMPsSSCP, adj.P.Val < 0.05)
dim(DMPsSSCPSigni)
#72 sondas significativas

UCSC<- strsplit(as.character(DMPsSSCPSigni$UCSC_RefGene_Name), ";")

UCSC <-unlist(UCSC)

Gencode<- strsplit(as.character(DMPsSSCPSigni$GencodeBasicV12_NAME), ";")

Gencode <-unlist(Gencode)

genes <- c(UCSC,Gencode)

genes <- unique(genes)

write.table(genes, file="genes.csv", sep=",", row.names=FALSE)





#Control vs proliferative DAD
#proliferative DAD: 2,6,7,19,21

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")
bVals <- readRDS("bVals.rds")
colnames(bVals) <- targets$EPIC.Sample.Name

MuestrasCP <- targets$EPIC.Sample.Name[targets$TISSUE %in% "CONTROL"] 
MuestrasCP <- c(MuestrasCP,"L2-P-Lu","L6-P-Lu","L7-P-Lu","L19-P-Lu","L21-P-Lu")

targetsCP <- targets[targets$EPIC.Sample.Name %in% c("L11-C-Lu","L13-C-Lu","L14-C-Lu","L17-C-Lu","L3-C-Lu","L4-C-Lu","L1-C-Lu","L15-C-Lu",
                                                     "L5-C-Lu","L7-C-Lu","L8-C-Lu","L10-C-Lu","L18-C-Lu","L24-C-Lu","L25-C-Lu","L19-C-Lu",
                                                     "L22-C-Lu","L23-C-Lu","L2-P-Lu","L6-P-Lu","L7-P-Lu","L19-P-Lu","L21-P-Lu"),]
bValsCP <- bVals[,colnames(bVals) %in% MuestrasCP]


# this is the factor of interest
tissue <- factor(targetsCP$TISSUE)
# this is the individual effect that we need to account for 
#gender <- factor(targets$GENDER) 
# use the above to create a design matrix
design <- model.matrix(~0+tissue, data=targetsCP) 
colnames(design) <- levels(tissue)
# fit the linear model 
#dim(mVals)
library(limma)
fit <- lmFit(bValsCP, design)
# create a contrast matrix for specific comparisons                     
contMatrix <- makeContrasts(CONTROL-COVID,
                            levels=design) 
contMatrix




fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))


#Este match es muy imprtante para que anote en orden
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850kSub <- ann850k[match(rownames(bVals),ann850k$Name),]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann850kSub)
DMPsSigni <- DMPs[DMPs$adj.P.Val < 0.05,]
dim(DMPsSigni)

bValsCP <- bVals[,colnames(bVals) %in% MuestrasCP]

colnames(bValsCP)

suma_renglones <- function (x) {
  y=integer(nrow(x))
  z=integer(nrow(x))
  r=integer(nrow(x))
  names <- colnames(x)
  for(i in 1:nrow(x)){
    for (j in 1:ncol(x)) {
      if(is.numeric(x[i,j]) && (names[j]=="L1-C-Lu"
                                || names[j]=="L3-C-Lu"
                                || names[j]=="L4-C-Lu"
                                || names[j]=="L5-C-Lu"
                                || names[j]=="L7-C-Lu"
                                || names[j]=="L8-C-Lu"
                                || names[j]=="L10-C-Lu"
                                || names[j]=="L11-C-Lu"
                                || names[j]=="L13-C-Lu"
                                || names[j]=="L14-C-Lu"
                                || names[j]=="L15-C-Lu"
                                || names[j]=="L17-C-Lu"
                                || names[j]=="L18-C-Lu"
                                || names[j]=="L19-C-Lu"
                                || names[j]=="L22-C-Lu"
                                || names[j]=="L23-C-Lu"
                                || names[j]=="L24-C-Lu"
                                || names[j]=="L25-C-Lu" )) 
      {
        y[i] = y[i]+(x[i,j])
      }else if (is.numeric(x[i,j]) && (names[j]=="L2-P-Lu"
                                       || names[j]=="L6-P-Lu"
                                       || names[j]=="L7-P-Lu"
                                       || names[j]=="L19-P-Lu"
                                       || names[j]=="L21-P-Lu")) {
        z[i] = z[i]+(x[i,j])
      }
      
    }
    dim(y)
    dim(x)
    y[i] <-y[i]/18
    z[i] <-z[i]/5
    r[i] = y[i]-z[i]
  }
  
  x = as.data.frame(x)
  x[, "mediaControl"] <- y
  x[, "mediaCovid"] <- z
  x[, "Diferencia"] <- r
  
  return(x)
}



BetaDiferenciaMedia = suma_renglones(bValsCP)

BetaDiferenciaMediaSS <- subset(BetaDiferenciaMedia, Diferencia > 0.2 | Diferencia < -0.2)
dim(BetaDiferenciaMediaSS)
#Me da 3652 sondas con diferencia entre medias mayor de 0.15
saveRDS(BetaDiferenciaMediaSS, "BetaDiferenciaMediaSSCP.rds")

#DMPsConVsAcuteDAD <- DMPs0.15Signi[rownames(DMPs0.15Signi) %in% rownames(BetaDiferenciaMediaSS),]

#BetaDiferenciaMediaSS <- readRDS("BetaDiferenciaMediaSSCP.rds")


#DMPsSSCP <- DMPs[match(rownames(BetaDiferenciaMediaSS),rownames(DMPs)),]
#(DMPsSSCP)

#DMPsSSCPSigni <- subset(DMPsSSCP, adj.P.Val < 0.05)
#dim(DMPsSSCPSigni)

DMPsSSCP <- DMPsSigni[match(rownames(BetaDiferenciaMediaSS),rownames(DMPsSigni)),]

#DMPsSSCP <- DMPsSigni[match(rownames(DMPsSigni),rownames(BetaDiferenciaMediaSS)),]

DMPsSSCP <-DMPsSSCP[!is.na(DMPsSSCP$chr),]
dim(DMPsSSCP)

BetaDiferenciaMediaSSCP <- BetaDiferenciaMediaSS[match(rownames(DMPsSSCP),rownames(BetaDiferenciaMediaSS)),]
write.table(BetaDiferenciaMediaSSCP, file="BetaDiferenciaMediaSSCP.csv", sep=",", row.names=TRUE)


UCSC<- strsplit(as.character(DMPsSSCP$UCSC_RefGene_Name), ";")

UCSC <-unlist(UCSC)

Gencode<- strsplit(as.character(DMPsSSCP$GencodeBasicV12_NAME), ";")

Gencode <-unlist(Gencode)

genes <- c(UCSC,Gencode)

genes <- unique(genes)

write.table(genes, file="genes.csv", sep=",", row.names=FALSE)



#Lista 10 genes
listaGenes <- topTable(fit2, num=10, coef=1, genelist=ann850kSub)
DMPs100 <- topTable(fit2, num=100, coef=1,sort.by="logFC", genelist=ann850kSub)

table(DMPs0.15Signi)


colnames(DMPsSSCP)
#tablaSondas
tablaSondas <- DMPsSSCP[,c("Name","chr","pos","GencodeBasicV12_NAME","GencodeBasicV12_Group","UCSC_RefGene_Name","UCSC_RefGene_Group","P.Value")] #Añadir columan posicion y cromosoma
write.table(tablaSondas, file="tablaSondas.csv", sep=",", row.names=TRUE)



#
DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
tablaSondas <- DMPs0.15Signi[,c("Name","chr","pos","GencodeBasicV12_NAME","GencodeBasicV12_Group","UCSC_RefGene_Name","UCSC_RefGene_Group","P.Value")] #Añadir columan posicion y cromosoma
write.table(tablaSondas, file="tablaSondas.csv", sep=",", row.names=TRUE)



UCSC<- strsplit(as.character(DMPs0.15Signi$GencodeBasicV12_Group), ";")

UCSC <- unique(UCSC)

UCSC <-unlist(UCSC)

Gencode<- strsplit(as.character(DMPs0.15Signi$UCSC_RefGene_Group), ";")

Gencode <-unlist(Gencode)

genes <- c(UCSC,Gencode)

genes <- unique(genes)

write.table(genes, file="genes.csv", sep=",", row.names=FALSE)



#GSE164013

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

library(tidyverse) #read_tsv
GSE164013 <- read_tsv("/Users/jparra/Documents/ProyectoCovid/data/GSE164013/Lung_Q3Norm_TargetCountMatrix.tsv")



#
BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")

BetaDiferenciaMediaSS <- BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% rownames(DMPs0.15Signi),]
write.table(BetaDiferenciaMediaSS, file="BetaDiferenciaMediaSS.csv", sep=",", row.names=FALSE)


#
getwd()
tablasondas84 <- read.csv("tablasondas84.csv")


#metilacion expresion graph para  PIP5K1B


library(tidyverse) #read_tsv
library(data.table) #as.data.table
library(ggpubr)  #ggplot #stat_cor

CosmicCompleteCL <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv")

betas_result450K<- read.csv("betas_result_450KAll1433.csv", skip=7) #Todo 450K

betas_resultEPIC <- read.csv("betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
betas_resultEPICSS <- betas_resultEPIC[colnames(betas_resultEPIC) %in% colnames(betas_result450K)]

betasAll <- rbind(betas_result450K,betas_resultEPICSS)
CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]
#k <- merge(CosmicSSAll,betasAll , by.x = "SAMPLE_NAME", by.y = "Sample.Name")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
DMPs0.15Signi$UCSC_RefGene_Name <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Nam, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_NAME <- sapply(strsplit(DMPs0.15Signi$GencodeBasicV12_NAME, ';', fixed=T), function(x) x[1])

DMPs0.15Signi$UCSC_RefGene_Group <- sapply(strsplit(DMPs0.15Signi$UCSC_RefGene_Group, ';', fixed=T), function(x) x[1])
DMPs0.15Signi$GencodeBasicV12_Group <- sapply(strsplit(DMPs0.15Signi$GencodeBasicV12_Group, ';', fixed=T), function(x) x[1])

#CosmicSS450K <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_result450K$Sample.Name,]
#CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% joinGenes,]

joinGenes <- read.csv("joinGenes.csv")
joinGenes <- joinGenes$x

#genesOcho <- c("IRF2BPL","NLRP3","ZNF608","EBF1","FSTL1","PIP5K1B","RAB11FIP5","ZBTB38")

#CosmicSS450K <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_result450K$Sample.Name,]
#CosmicSS450KGenes <- CosmicSS450K[CosmicSS450K$GENE_NAME %in% joinGenes,]

#filtramos lineas celulares
CosmicSSEPIC <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name,]
#Despues filtramos por genes. En mi caso obtengo 10360 lineas.
CosmicSSEPICGenes <- CosmicSSEPIC[CosmicSSEPIC$GENE_NAME %in% joinGenes,]

#Tambien podriamos filtrar directamente por genes pero obtenog muchas lienas celulares que no voy a usar. Tengo 71780.
#CosmicSSEPICGenes <- CosmicCompleteCL[CosmicCompleteCL$GENE_NAME %in% joinGenes,]


 
#betasAll[,colnames(betasAll) %in% "cg26818401"]
#betas_resultEPIC[,colnames(betas_resultEPIC) %in% "cg26818401"]
#betas_result450K[,colnames(betas_result450K) %in% "cg26818401"]


pdf("graficas expresion-metilacion/GraficaEpicAllGencode0407.pdf")
#pdf("graficas expresion-metilacion/GraficaEpic+450KAllGencode0407.pdf")

#library(ggpubr)  #ggplot #stat_cor
#library(ggrepel)  #geom_text_repel
i <- "UNC80" 

for (i in joinGenes) {
  #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
  if(i %in% CosmicSSEPIC$GENE_NAME){
    
    rm(sonda)
    
    j <- assign(paste0("gen", i), CosmicSSEPICGenes[CosmicSSEPICGenes$GENE_NAME %in% i,])
    #j <- assign(paste0("gen", i), CosmicSS450KGenes[CosmicSS450KGenes$GENE_NAME %in% i,])
    #j <- assign(paste0("gen", i), CosmicSSAll[CosmicSSAll$GENE_NAME %in% i,])

    sonda <- DMPs0.15Signi[DMPs0.15Signi$UCSC_RefGene_Name %in% i, "Name" ]
    #sonda <- DMPs0.15Signi[DMPs0.15Signi$GencodeBasicV12_NAME %in% i, "Name" ]
    
    sonda <- sonda[1]
    
    if(!is_empty(sonda) && sonda  %in% colnames(betas_resultEPIC)){
      expresionGen <- betas_resultEPIC[betas_resultEPIC$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      # expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
      expresionGen <- as.data.table(expresionGen)
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
      #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]
      
      plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
        labs(x=sonda) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        #theme(panel.background = element_rect(fill = "#67c9ff")) +
        theme_bw() +
        #geom_text_repel(label=k$SAMPLE_NAME) +
        #geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
        ggtitle(paste("Expression-methylation graph for the gene:",i,"that is found in:",variablePromotor))
    
      print(plot)
      
    }
     }
   } #a
  
# }

dev.off()










#TCGA

#VER sondas de IRF2BPL para encontrar ekivalene
#Hacxer emtilacione xpresionc on TCGA.


ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

ann850kIRF <- ann850k[ann850k$GencodeBasicV12_NAME %in% "IRF2BPL",]

ann850kIRF <- as.data.frame(ann850kIRF)

ann850k14 <- ann850k[ann850k$chr %in% "chr14",]

ann850k14 <- as.data.frame(ann850k14)


ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k14 <- ann450k[ann450k$chr %in% "chr14",]
ann450k14 <- as.data.frame(ann450k14)

ann450kIRF <- ann450k[ann450k$UCSC_RefGene_Name %in% "IRF2BPL",]




BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
#hg38 no funciona
install.packages("IlluminaHumanMethylationEPICanno.ilm10b5.hg38")

library(devtools)
install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#esta anotacion tampoco coincide con TCGA

ann850kIRF <- ann850k[ann850k$UCSC_RefGene_Name %in% "IRF2BPL",]
ann850kIRF <- as.data.frame(ann850kIRF)


install.packages("ghit")
library(ghit)
install.packages("tibble")
library(ghit)





vectorGenesEnrichr <- c("PLEKHG1","PLXND1","CTNND1","EBF1","ATP11A","THADA","IRF2BPL","NCOR2","SKI","FMNL3","ZNF608","CDK6","ADORA3","RAPGEF2",
                        "RAB11FIP5","PLEKHG1","ADCY9","CDK6","MAST4","CTNND1","CD38","LDB3","NREP","TRAK1","FSTL1","NCOR2","SKI","SPON2","ZNF608","ADCY9","LDB3","ATP11A",
                        "NREP","GFRA2","DSTN","TCF4","NUPR1","NREP","FSTL1","PLXND1","CTNND1","RAPGEF2","BIRC6","ATP11A","TRAK1","SPTBN1","PIAS1","PLEKHG1","TPPP3","ADCY9",
                        "MAST4","MCF2L","RAPGEF2","TCF4","SPTBN1","DSTN","TCF4","NUPR1","NREP","FSTL1","KIAA0930","NCOR2","CTNND1","RAPGEF2","SPTBN1","IRF2BPL",
                        "ZNF608","CDK6","ZBTB38","RAPGEF2","BIRC6","ATP11A","TRAK1","PIAS1","ZNF608","CDK6","HOXA3","ITPR1","PIP5K1B","TCF4","TRAM2","IRF2BPL","FMNL3",
                        "PLEKHG1","ZNF608","PLXND1","GPX3","CTNND1","RAPGEF2","ADCY9","MAST4","MCF2L","HOXA3","CD38","RAPGEF2","TCF4") 

vectorGenesEnrichrUnique <- unique(vectorGenesEnrichr)
vectorGenesEnrichrUnique <- as.data.frame(vectorGenesEnrichrUnique)
write.table(vectorGenesEnrichrUnique, file="vectorGenesEnrichrUnique.csv", sep=",", row.names=FALSE)



count(vectorGenesEnrichrUnique, vectorGenesEnrichrUnique[1])

vectorGenesEnrichrUnique <- vectorGenesEnrichrUnique[[1]]

library(dplyr)

count(vectorGenesEnrichrUnique, vectorGenesEnrichrUnique[1])


count(vectorGenesEnrichrUnique, "FSTL1")

library("stringr")                       # Load stringr package
stringr::str_count(vectorGenesEnrichrUnique, "FSTL1")


vectorGenesEnrichr <- as.data.frame(vectorGenesEnrichr)
write.table(vectorGenesEnrichr, file="vectorGenesEnrichr.csv", sep=",", row.names=FALSE)



#heatmap 119
library(heatmap3)
setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

bVals2 <- readRDS("bVals.rds")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")
DMPs0.15Signi <- DMPs0.15Signi120

bVals0.15Signi <- bVals[match(rownames(DMPs0.15Signi),rownames(bVals)),]
colnames(bVals0.15Signi) <- targets$EPIC.Sample.Name
heatmap3(bVals0.15Signi, ColSideCut=10)
heatmap3(bVals0.15Signi)



DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

bValsSS0.15 <- bVals[rownames(bVals) %in% rownames(DMPs0.15Signi),]

bVals0.15Signi <- bVals[match(rownames(DMPs0.15Signi),rownames(bVals)),]
colnames(bVals0.15Signi) <- targets$EPIC.Sample.Name

library(sjmisc) #str_contains
x <- c()
for (i in colnames(bVals0.15Signi)){
  if( str_contains(i, "P")){
    x <- c(x,"blue")
  }else{
    x <- c(x,"yellow")
  }
  
}


#heatmap3(bVals0.15Signi, ColSideCut=10, ColSideColors=x, method = "ward.D2") 


heatmap3(bVals0.15Signi, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 


ann850kSS <- ann850k[,c("Name","chr","Islands_Name","Relation_to_Island","UCSC_RefGene_Name", "UCSC_RefGene_Group", "GencodeBasicV12_NAME","GencodeBasicV12_Group" )]
DMPs0.15SigniSS <- DMPs0.15Signi[,c("Name","chr","Islands_Name","Relation_to_Island","UCSC_RefGene_Name", "UCSC_RefGene_Group", "GencodeBasicV12_NAME","GencodeBasicV12_Group" )]
write.table(DMPs0.15SigniSS, file="DMPs0.15SigniSS.csv", sep=",", row.names=FALSE)




#tabla sondas por gen:
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

nombreGenesUSDCSplit<- strsplit(as.character(ann850k$UCSC_RefGene_Name), ";")

#nombreGenesUSDCSplitUnicos <- unique(nombreGenesUSDCSplit)

nombreGenesUSDCSplitUnicosUnlist <-unlist(nombreGenesUSDCSplit)

nombreGenesUSDCSplitUnicos <- unique(nombreGenesUSDCSplitUnicosUnlist)

length(nombreGenesUSDCSplitUnicos)
#27364


nombreGenesUSDCSplitUnicos[1]

#sondas <- ann850k$Name[ann850k$UCSC_RefGene_Name %in% "MID1IP1"]

GenesYSondas <-data.frame()
count <-1
#gen <- "MID1IP1"
for (gen in nombreGenesUSDCSplitUnicos){
  sondas <- ann850k$Name[grepl(gen, ann850k$UCSC_RefGene_Name)]
  sondas<- as.character(sondas) 
  GenesYSondas[count,1] <- gen
  GenesYSondas[count,2] <- concat(sondas,sep=";")
  count <- count+1
  
}




sondas <- ann850k$Name[grepl("MID1IP1", ann850k$UCSC_RefGene_Name)]
sondas<- as.character(sondas) 
GenesYSondas[1,1] <- nombreGenesUSDCSplitUnicos[1]
GenesYSondas[1,2] <- concat(sondas,sep=";")

write.table(GenesYSondas, file="GenesYSondas.csv", sep=",", row.names=FALSE)




sondasCon<- c(sondas)
sondasCon[1]<-sapply(sondas, paste, collapse = ";")

install.packages("bazar")
library(bazar)

sondasCon<- concat(sondas,sep=";")
ann850k$UCSC_RefGene_Name[ann850k$Name %in% "cg13298827"]
ann850kDF <- as.data.frame(ann850k)

ann850k$UCSC_RefGene_Name

ann850kSondasYGenes <- ann850k[,c("Name","UCSC_RefGene_Name")]
ann850kSondasYGenes <- as.data.frame(ann850kSondasYGenes)
write.table(ann850kSondasYGenes, file="ann850kSondasYGenes.csv", sep=",", row.names=FALSE)


#heatmap para linux

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

library(heatmap3)

bValsRSS <- readRDS("bValsSS.rds")


#library(sjmisc) #str_contains
#x <- c()
#for (i in colnames(bValsRSS)){
#  if( str_contains(i, "P")){
#    x <- c(x,"blue")
#  }else{
#    x <- c(x,"yellow")
#  }
#}



heatmap3(bValsRSS, ColSideCut=1, ColSideColors=x, method = "ward.D2") 

heatmap3(bValsRSS, method = "ward.D2")

ann850k$UCSC_RefGene_Group



#RAPGEF2 

sondas <- ann850k$Name[grepl("RAPGEF2", ann850k$GencodeBasicV12_NAME)]

ann850kRAPGEF2 <- ann850k[ann850k$Name %in% sondas,c("pos","GencodeBasicV12_NAME","GencodeBasicV12_Group")]
ann850kRAPGEF2 <- as.data.frame(ann850kRAPGEF2)


BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")
BetaDiferenciaMedia$Name <- rownames(BetaDiferenciaMedia)
BetaDiferenciaMediaRAPGEF2  <- BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% sondas,c("Diferencia","Name")]
ann850kRAPGEF2$Diferencia <- BetaDiferenciaMediaRAPGEF2[match(rownames(ann850kRAPGEF2),rownames(BetaDiferenciaMediaRAPGEF2)),]



#UCSC

sondas <- ann850k$Name[grepl("NLRP3", ann850k$UCSC_RefGene_Name)]

ann850kRAPGEF2 <- ann850k[ann850k$Name %in% sondas,c("pos","UCSC_RefGene_Name","UCSC_RefGene_Group")]
ann850kRAPGEF2 <- as.data.frame(ann850kRAPGEF2)


BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")
BetaDiferenciaMedia$Name <- rownames(BetaDiferenciaMedia)
BetaDiferenciaMediaRAPGEF2  <- BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% sondas,c("Diferencia","Name")]
ann850kRAPGEF2$Diferencia <- BetaDiferenciaMediaRAPGEF2[match(rownames(ann850kRAPGEF2),rownames(BetaDiferenciaMediaRAPGEF2)),]


DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")



#NLRP3


sondas <- ann850k$Name[grepl("NLRP3", ann850k$GencodeBasicV12_NAME)]

ann850kNLRP3 <- ann850k[ann850k$Name %in% sondas,c("pos","GencodeBasicV12_NAME","GencodeBasicV12_Group")]
ann850kNLRP3 <- as.data.frame(ann850kNLRP3)


BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")
BetaDiferenciaMedia$Name <- rownames(BetaDiferenciaMedia)
BetaDiferenciaMediaNLRP3  <- BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% sondas,c("Diferencia","Name")]
ann850kNLRP3$Diferencia <- BetaDiferenciaMediaNLRP3[match(rownames(ann850kNLRP3),rownames(BetaDiferenciaMediaNLRP3)),]


DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

a <-rownames(ann850k)
a <- as.data.frame(a)
#UCSC

sondas <- ann850k$Name[grepl("NLRP3", ann850k$UCSC_RefGene_Name)]

ann850kNLRP3 <- ann850k[ann850k$Name %in% sondas,c("pos","UCSC_RefGene_Name","UCSC_RefGene_Group")]
ann850kNLRP3 <- as.data.frame(ann850kNLRP3)


BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")
BetaDiferenciaMedia$Name <- rownames(BetaDiferenciaMedia)
BetaDiferenciaMediaNLRP3  <- BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% sondas,c("Diferencia","Name")]
ann850kNLRP3$Diferencia <- BetaDiferenciaMediaNLRP3[match(rownames(ann850kNLRP3),rownames(ann850kNLRP3)),]


DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")



ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)




#CrossReactive:

remotes::install_github("markgene/maxprobes")

library(maxprobes)
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")
length(xloci)
xloci<- unlist(xloci)
xloci<- unique(xloci)
xloci <- as.data.frame(xloci)
write.table(xloci, file="crossReactiveEpic.csv", sep=",")



xloci <- maxprobes::xreactive_probes(array_type = "450K")
length(xloci)
xloci<- unlist(xloci)
xloci<- unique(xloci)
xloci <- as.data.frame(xloci)
write.table(xloci, file="crossReactive450K.csv", sep=",")


rgSet <- rgSet[minfi::featureNames(rgSet)[!(minfi::featureNames(rgSet) %in% bad_pos)], ]



ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


#crossReactive
library(maxprobes)
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")
#xloci <- maxprobes::xreactive_probes(array_type = "450K")

xloci<- unlist(xloci)
xloci<- unique(xloci)

ann850kSS <- ann850k[!ann850k$Name %in%  xloci,]


#Para filtrar snps
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)


#Para filtrar ChrX-Y
ann850kSS <- ann850k[!ann850k$chr %in%   c("chrX","chrY"),]





#120 o 119

getwd()

setwd("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias")

BetaDiferenciaMedia <- readRDS("BetaDiferenciaMedia.rds")

DMPs <- readRDS( "DMPs.rds")


#BetaDiferenciaMedia0.15 <- subset(BetaDiferenciaMedia, Diferencia > 0.15 | Diferencia < -0.15)
BetaDiferenciaMedia0.15 <- BetaDiferenciaMedia[BetaDiferenciaMedia$Diferencia > 0.15 | BetaDiferenciaMedia$Diferencia < -0.15,]



#BetaDiferenciaMedia[rownames(BetaDiferenciaMedia) %in% "cg27496650", ]
#BetaDiferenciaMedia0.15[rownames(BetaDiferenciaMedia0.15) %in% "cg27496650", ]

DMPs0.15 <- DMPs[match(rownames(BetaDiferenciaMedia0.15),rownames(DMPs)),]


#DMPs0.15Signi <- subset(DMPs0.15, adj.P.Val < 0.05) #Los p.valores son difentees antes y despues del filtrado!!!! no usar subset!


DMPs[DMPs$Name %in% "cg27496650", ]


DMPs0.15[DMPs0.15$Name %in% "cg27496650", ]
DMPs0.15Signi[DMPs0.15Signi$Name %in% "cg27496650", ]

DMPs0.15[DMPs0.15$Name %in% "cg08884752", ]
DMPs0.15Signi[DMPs0.15Signi$Name %in% "cg08884752", ]
DMPs0.15Signi120[DMPs0.15Signi120$Name %in% "cg08884752", ]



DMPs0.15Signi <- DMPs0.15[DMPs0.15$adj.P.Val < 0.05,]
DMPs0.15Signi120 <- DMPs0.15Signi


DMPs0.15Signi120[!DMPs0.15Signi120$Name %in% DMPs0.15Signi$Name,]

DMPs0.15Signi120 <- readRDS("DMPs0.15Signi120.rds")

DMPs0.15Signi120[!DMPs0.15Signi120$Name %in% DMPs0.15Signi$Name,]

#Graficas met-expr









#Tabla2



BetaDiferenciaMedia <- read.csv("BetaDiferenciaMedia.csv")

DMPs0.15Signi120[DMPs0.15Signi120$Name %in% "cg27496650", ]
#write.table(BetaDiferenciaMedia, file="BetaDiferenciaMedia.csv", sep=",", row.names=TRUE)

DMPs0.15Signi120 <- readRDS("DMPs0.15Signi120.rds")

BetaDiferenciaMedia[BetaDiferenciaMedia$X %in% "cg27496650", ]






BetaDiferenciaMedia <- read.csv("BetaDiferenciaMedia.csv")
BetaDiferenciaMedia <- read.csv("BetaDiferenciaMedia1.csv")


sondasYnombres <- DMPs0.15Signi120[,c("GencodeBasicV12_NAME","UCSC_RefGene_Name","Name","chr","pos","adj.P.Val")] #Añadir columan posicion y cromosoma


BetaDiferenciaMedia120 <- BetaDiferenciaMedia[ rownames(BetaDiferenciaMedia) %in%  DMPs0.15Signi120$Name,  ]

DMPs0.15Signi120 <- DMPs0.15Signi120[match((DMPs0.15Signi120$Name),rownames(BetaDiferenciaMedia120)),] 

sondasYnombres <- c(DMPs0.15Signi120,BetaDiferenciaMedia120)
#sondasYnombres$GencodeCompV12_NAME
sondasYnombres <- as.data.frame(sondasYnombres)

sondasYnombres <-sondasYnombres[,c("Name","chr","pos","GencodeBasicV12_NAME", "GencodeBasicV12_Group",  "adj.P.Val","mediaControl","mediaCovid","Diferencia")]
write.table(sondasYnombres, file="Tabla2Gencode.csv", sep=",", row.names=TRUE)
sondasYnombres$GencodeBasicV12_Group

#gometh https://rdrr.io/bioc/missMethyl/man/gometh.html

BiocManager::install("missMethyl")

library(missMethyl)

sigcpgs <- DMPs0.15Signi$Name
#sigcpgs <- as.array(sigcpgs)
allcpgs <- rownames(bVals)
#allcpgs <- as.array(allcpgs)

gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "KEGG", 
              array.type =  "EPIC")
help(gometh)


#Probar con GSEA o inmunne gene set

result <- gsameth(sig.cpg = sigcpgs, all.cpg = allcpgs, 
                  collection = go$idList, array.type = array.type, 
                  plot.bias = plot.bias, prior.prob = prior.prob, 
                  anno = anno, equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                  genomic.features = genomic.features, sig.genes = sig.genes)


result <- gsameth(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO")


listaGenesGencode <- DMPs0.15Signi$GencodeBasicV12_NAME
listaGenesGencode <- na.omit(listaGenesGencode)
listaGenesGencode<- strsplit(as.character(listaGenesGencode), ";")
listaGenesGencode <-unlist(listaGenesGencode)
listaGenesGencode <- unique(listaGenesGencode)
write.table(listaGenesGencode, file="listaGenesGencode.csv", sep=",", row.names=TRUE)

bVals0.15Signi <- as.data.frame(bVals0.15Signi)
a <-bVals0.15Signi[rownames(bVals0.15Signi) %in% "cg02471028",]

bVals <- as.data.frame(bVals)
b <-bVals[rownames(bVals) %in% "cg03719830",]


library(ComplexHeatmap)
library(Cairo)
palette_greenred <- colorRampPalette( c("green", "black", "red"), space="rgb")(75)

Heatmap(as.matrix(bVals0.15Signi), name = "mat",cluster_rows=T,use_raster=TRUE,
        show_row_names=FALSE,show_column_names=FALSE,
        raster_device='CairoPNG', col = palette_greenred,
           row_dend_reorder = TRUE,show_row_dend = F)




saveRDS(DMPs0.15Signi, "DMPs0.15Signi.rds")

DMPs0.15Signi <- readRDS("DMPs0.15Signi120.rds")


remove.packages("readr")
install.packages("readr")



#bVals0.15Signi
bVals <- readRDS("bVals.rds")

DMPs0.15Signi <- readRDS("DMPs0.15Signi.rds")

bVals0.15Signi <- bVals[match(rownames(DMPs0.15Signi),rownames(bVals)),]



####Diferenctes resultados por el preprocesado.
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
install.packages("readr")
library(readr)
library(minfi)

dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")

#rgSet <- read.metharray.exp(base="/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/IdatsAutopsiasFinish34",  targets=targets)
rgSet <- readRDS("rgSet.rds")



#minfi
detP <- detectionP(rgSet)




#Filtro MUESTRAS con media de p-valor menor de 0.05.
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

targets <- targets[keep,]
detP <- detP[,keep]
dim(detP)


# normalizacion 

#Como lo usan en el script:     mSetSq <- minfi::preprocessNoob(rgSet, dyeMethod = "single")

#mSetSq <- preprocessQuantile(rgSet) 
##Resultado:  865859 sondas

mSetSq <- preprocessNoob(rgSet) 

mSetSq <- preprocessNoob(rgSet, dyeMethod = "single") 

mSetSq <- preprocessFunnorm(rgSet) 


#preprocessRaw, preprocessFunnorm ,preprocessNoob



# ensure probes are in the same order in the mSetSq and detP objects    
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 


#Irene:Quitamos las SONDAS que tengan una p.valor media mayor de 0.01
bad_pos <- row.names(as.data.frame(detP))[rowMeans(detP) > 0.01]
rgSet <- rgSet[minfi::featureNames(rgSet)[!(minfi::featureNames(rgSet) %in% bad_pos)], ]

mSetSq <- mSetSq[minfi::featureNames(mSetSq)[!(minfi::featureNames(mSetSq) %in% bad_pos)], ]
#Resultado:  823254 sondas


#Filtrado cromosomas sexuales
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

keep <- !(featureNames(mSetSq) %in% ann850k$Name[ann850k$chr %in%   c("chrX","chrY")])

table(keep)
mSetSqFlt <- mSetSq[keep,]
##Resultado:  806600 sondas



# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
#Resultado:  779078 sondas



#saveRDS(mSetSqFlt, "mSetSqFlt.rds") 
#mSetSqFlt <- readRDS("mSetSqFlt.rds")





#Filtramos sondas crossReactives
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")

xloci<- unlist(xloci)
xloci<- unique(xloci)

keep <- !(featureNames(mSetSqFlt) %in% xloci)
table(keep)

mSetSqFlt <- mSetSqFlt[keep,]
#Resultado:  739071 sondas


# calculate beta for statistical analysis

bVals <- getBeta(mSetSqFlt)
bValsFunNorm <- as.data.frame(bVals)
bValsQuantileEjemplo1 <- bValsQuantile[rownames(bValsQuantile) %in% "cg11189357",colnames(bValsQuantile) %in% "204776850059_R08C01"] 
bValsNoobEjemploS1 <- bValsNoobSimple[rownames(bValsNoobSimple) %in% "cg11189357",colnames(bValsNoobSimple) %in% "204776850059_R08C01"] 
bValsFunNormEjemplo1 <- bValsFunNorm[rownames(bValsFunNorm) %in% "cg11189357",colnames(bValsFunNorm) %in% "204776850059_R08C01"] 


bValsQuantileEjemplo2 <- bValsQuantile[rownames(bValsQuantile) %in% "cg02471028",colnames(bValsQuantile) %in% "204860850155_R03C01"] 
bValsNoobEjemploS2 <- bValsNoobSimple[rownames(bValsNoobSimple) %in% "cg02471028",colnames(bValsNoobSimple) %in% "204860850155_R03C01"]
#bValsNoobEjemploS1 <- bValsNoobSimple[rownames(bValsNoobSimple) %in% "cg02471028",colnames(bValsNoobSimple) %in% "204860850155_R03C01"] 
bValsFunNormEjemplo2 <- bValsFunNorm[rownames(bValsFunNorm) %in% "cg02471028",colnames(bValsFunNorm) %in% "204860850155_R03C01"] 

bValsQuantileEjemplo3 <- bValsQuantile[rownames(bValsQuantile) %in% "cg02471028",colnames(bValsQuantile) %in% "204776850059_R01C01"] 


detP[rownames(detP) %in% "cg02471028",colnames(detP) %in% "204860850155_R03C01"] 





####SCRIPT

#probe_detP_max <- 0.01

dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script"
dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"

targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")

targets[1, ]
#process_sample_methylarray(rgSet,targets)
#process_sample_methylarray(1:34,targets,"/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script",0.01 )
process_sample_methylarray(1:34,targets,"/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script",0.25,0.01 )

sample1 <- read.csv(file="/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script/x204776850059_R01C01.csv")

scriptEjemplo3 <- sample1[sample1$TargetID %in% "cg02471028","x204860850155_R01C01"] ##Casi el mismo valor que en methylationDB.


  for (val in 1:34) {
    process_sample_methylarray(val,targets,"/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script2",0.25,0.01 )
  }



sample2 <- read.csv(file="/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script2/x204860850155_R01C01.csv")

scriptEjemplo3.2 <- sample2[sample2$TargetID %in% "cg02471028","x204860850155_R01C01"] ##Casi el mismo valor que en methylationDB.




bValsNoobEjemploS2 <- bValsNoobSimple[rownames(bValsNoobSimple) %in% "cg08884752",]

                           

#Ejemplo linea celular.

dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script3/"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop.csv")
targets[1, ]
targets$Basename <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script3/idats/204088090136_R06C01"

process_sample_methylarray(1,targets,"/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/script3/idats",0.25,0.01 )



#testeo
bVals[rownames(bVals) %in% "cg02471028","204776850059_R07C01"] ##Casi el mismo valor que en methylationDB.





#Graficas metilacion expresion el pulmon y hematopoyeticas.


#NLRP3. sonda: 18126557

#pulmon

betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
betas450KLung <- betas_result450K[betas_result450K$Site.Primary %in% "Lung",] #223 muestras

betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
betasEPICLung <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Lung",] #223 muestras
betasEPICLungSS <- betasEPICLung[colnames(betasEPICLung) %in% colnames(betas450KLung)]
#betas_resultEPICSS tiene 373 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.

#NLRP3. sonda: 18126557
colnames(betasEPICLungSS) %in% "cg18126557"

betasAll <- rbind(betas450KLung,betasEPICLungSS)
#244 muetras de pulmn entre 450 y epic

CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]



      
      j <- assign(paste0("gen", "NLRP3"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "NLRP3",])
      
       expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg18126557")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")

      plot <- ggplot(k,aes(x = k[ ,"cg18126557"],y = Z_SCORE)) +
        labs(x="cg18126557") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","NLRP3","that is found in:","3'UTR"))
      
      print(plot)
      


      
      
      #Hematopoyeticas
      unique(betas_result450K$Site.Primary)
      betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
      betas450KHema <- betas_result450K[betas_result450K$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      
      betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
      betasEPICHema <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      betasEPICHemaSS <- betasEPICHema[colnames(betasEPICHema) %in% colnames(betas450KHema)]
      #betas_resultEPICSS tiene 168 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.
      
      #NLRP3. sonda: 18126557
      colnames(betasEPICLungSS) %in% "cg18126557"
      
      betasAll <- rbind(betas450KHema,betasEPICHemaSS)
      #363 muetras hemato entre 450 y epic
      
      CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]
      
      
      
      
      j <- assign(paste0("gen", "NLRP3"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "NLRP3",])
      
      expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg18126557")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      plot <- ggplot(k,aes(x = k[ ,"cg18126557"],y = Z_SCORE)) +
        labs(x="cg18126557") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","NLRP3","that is found in:","3'UTR"))
      
      print(plot)
      
      

      
      #REVISAR SI ESTOS 4 GENES TIENEN MAS DE UNA SONDA, en las 120 casi seguro que no
      
      
      #ZNF608 sonda: cg06055229
      
      #pulmon
      
      betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
      betas450KLung <- betas_result450K[betas_result450K$Site.Primary %in% "Lung",] #223 muestras
      
      betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
      betasEPICLung <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Lung",] #223 muestras
      betasEPICLungSS <- betasEPICLung[colnames(betasEPICLung) %in% colnames(betas450KLung)]
      #betas_resultEPICSS tiene 373 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.
      
      #NLRP3. sonda: 18126557
      colnames(betasEPICLungSS) %in% "cg06055229"
      
      betasAll <- rbind(betas450KLung,betasEPICLungSS)
      #244 muetras de pulmn entre 450 y epic
      
      CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]
      
      
      
      
      j <- assign(paste0("gen", "ZNF608"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "ZNF608",])
      
      expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg06055229")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      plot <- ggplot(k,aes(x = k[ ,"cg06055229"],y = Z_SCORE)) +
        labs(x="cg06055229") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","ZNF608","that is found in:","5'UTR"))
      
      print(plot)
      
      
      
      
      
      #Hematopoyeticas
      unique(betas_result450K$Site.Primary)
      betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
      betas450KHema <- betas_result450K[betas_result450K$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      
      betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
      betasEPICHema <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      betasEPICHemaSS <- betasEPICHema[colnames(betasEPICHema) %in% colnames(betas450KHema)]
      #betas_resultEPICSS tiene 373 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.

      
      betasAll <- rbind(betas450KHema,betasEPICHemaSS)
      #363 muetras hemato entre 450 y epic
      
      CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]
      
      
      
      
      j <- assign(paste0("gen", "ZNF608"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "ZNF608",])
      
      expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg06055229")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      plot <- ggplot(k,aes(x = k[ ,"cg06055229"],y = Z_SCORE)) +
        labs(x="cg06055229") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","ZNF608","that is found in:","5'UTR"))
      
      print(plot)
      
      
      
      
      #IRF2BPL sonda: cg13493526
      
      #pulmon
      
      betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
      betas450KLung <- betas_result450K[betas_result450K$Site.Primary %in% "Lung",] #223 muestras
      
      betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
      betasEPICLung <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Lung",] #223 muestras
      betasEPICLungSS <- betasEPICLung[colnames(betasEPICLung) %in% colnames(betas450KLung)]
      #betas_resultEPICSS tiene 373 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.
      
      colnames(betasEPICLungSS) %in% "cg13493526"
      
      betasAll <- rbind(betas450KLung,betasEPICLungSS)
      #244 muetras de pulmn entre 450 y epic
      
      CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]
      
      
      
      
      j <- assign(paste0("gen", "IRF2BPL"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "IRF2BPL",])
      
      expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg13493526")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      plot <- ggplot(k,aes(x = k[ ,"cg13493526"],y = Z_SCORE)) +
        labs(x="cg13493526") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","IRF2BPL","that is found in:","1stExon"))
      
      print(plot)
      
      
      
      
      
      #Hematopoyeticas
      unique(betas_result450K$Site.Primary)
      betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
      betas450KHema <- betas_result450K[betas_result450K$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      
      betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
      betasEPICHema <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      betasEPICHemaSS <- betasEPICHema[colnames(betasEPICHema) %in% colnames(betas450KHema)]
      #betas_resultEPICSS tiene 373 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.
      
      
      betasAll <- rbind(betas450KHema,betasEPICHemaSS)
      #363 muetras hemato entre 450 y epic
      
      CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasAll$Sample.Name,]
      
      
      
      
      j <- assign(paste0("gen", "IRF2BPL"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "IRF2BPL",])
      
      expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg13493526")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      plot <- ggplot(k,aes(x = k[ ,"cg13493526"],y = Z_SCORE)) +
        labs(x="cg13493526") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","IRF2BPL","that is found in:","1stExon"))
      
      print(plot)
      
      
      
      
      #RAPGEF2 sonda: cg02717112
      
      #pulmon
      
      betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
      betas450KLung <- betas_result450K[betas_result450K$Site.Primary %in% "Lung",] #223 muestras
      
      betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
      betasEPICLung <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Lung",] #223 muestras
      betasEPICLungSS <- betasEPICLung[colnames(betasEPICLung) %in% colnames(betas450KLung)]
      #betas_resultEPICSS tiene 373 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.
      
      colnames(betasEPICLung) %in% "cg02717112"
      colnames(betas450KLung) %in% "cg02717112"
      #Esta sonda esta en epic pero no en 450K.
      
     # betasAll <- rbind(betas450KLung,betasEPICLungSS)
      #244 muetras de pulmn entre 450 y epic
      
      CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasEPICLung$Sample.Name,]
      
      
      
      
      j <- assign(paste0("gen", "RAPGEF2"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "RAPGEF2",])
      
      expresionGen <- betasEPICLung[betasEPICLung$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg02717112")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      plot <- ggplot(k,aes(x = k[ ,"cg02717112"],y = Z_SCORE)) +
        labs(x="cg02717112") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","RAPGEF2","that is found in:","TSS200"))
      
      print(plot)
      
      
      
      
      
      #Hematopoyeticas
      unique(betas_result450K$Site.Primary)
      betas_result450K<- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_result_450KAll1433.csv", skip=7) #Todo 450K
      betas450KHema <- betas_result450K[betas_result450K$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      
      betas_resultEPIC <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/betas_resultAllCellLinesEPIC373.csv", skip=7) #Todo EPIC
      betasEPICHema <- betas_resultEPIC[betas_resultEPIC$Site.Primary %in% "Haematopoietic and Lymphoid",] #223 muestras
      betasEPICHemaSS <- betasEPICHema[colnames(betasEPICHema) %in% colnames(betas450KHema)]
      #betas_resultEPICSS tiene 373 observaciones para 57 sondas,de las 120 solo 52 estan en 450K.
      
      
     # betasAll <- rbind(betas450KHema,betasEPICHemaSS)
      #363 muetras hemato entre 450 y epic
      
      CosmicSSAll <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betasEPICHema$Sample.Name,]
      
      
      
      j <- assign(paste0("gen", "RAPGEF2"), CosmicSSAll[CosmicSSAll$GENE_NAME %in% "RAPGEF2",])
      
      expresionGen <- betasEPICHema[betasEPICHema$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", "cg02717112")]
      expresionGen <- as.data.table(expresionGen)
      #Si tengo mas de un valor para una linea celular hago la media.
      expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
      
      k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
      
      plot <- ggplot(k,aes(x = k[ ,"cg02717112"],y = Z_SCORE)) +
        labs(x="cg02717112") +
        geom_point() +
        geom_smooth(method='lm') +
        stat_cor(method = "spearman" ) +
        theme_bw() +
        ggtitle(paste("Expression-methylation graph for the gene:","RAPGEF2","that is found in:","TSS200"))
      
      print(plot)
      
      
      
      
      #Control de calidad
      rgSet <- readRDS("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/rgSet.rds")
      dataDirectory <- "/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias"
      targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheetAutop34.csv")
      qcReport(rgSet, sampNames=targets$SAMPLE.ID, sampGroups=targets$TISSUE, 
               pdf="qcReport.pdf")
      
      
      
      
      #Resultados GSEA 73
      rm(ResultsGSEA)
      ResultsGSEA <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/resultGSEA73.csv", TRUE, sep = ",")
      colnames(ResultsGSEA)
      plot <-ggplot(ResultsGSEA, aes(x=Gene.Set.Name,  y = X..Genes.in.Gene.Set..K.))
      plot
      
      
      barplot_GSEA <- function(data = data, title = "") {
        data$Gene.Set.Name <- gsub("^.*?_","", data$Gene.Set.Name)
        data$Gene.Set.Name <- factor(data$Gene.Set.Name, levels = data$Gene.Set.Name[order(as.numeric(data$X..Genes.in.Gene.Set..K.), decreasing = FALSE)])
        ggplot(data, aes_string(x = "X..Genes.in.Overlap..k.", y = "Gene.Set.Name", fill = "p.value")) +
          geom_col()+
          theme_bw() +
          scale_fill_continuous(low="red", high="red4", name = "p.value",
                                guide=guide_colorbar(reverse=TRUE)) +
          theme(axis.text = element_text(size = 12)) +
          ggtitle(title) +
          xlab(NULL) + ylab(NULL)
      }
      
      library(ggplot2)
      barplot_GSEA(ResultsGSEA)
      plot
      
      ResultsGSEA$p.value

      
      
      
      
 ResultsGSEA <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/resultGSEA73.csv", TRUE, sep = ",")
 ResultsGSEA$p.value<- as.numeric(ResultsGSEA)
 
 barplot_GSEA <- function(data = data, title = "") {
   #ResultsGSEA$Gene.Set.Name <- gsub("^.*?_","", ResultsGSEA$Gene.Set.Name)
   ResultsGSEA$Gene.Set.Name <- factor(ResultsGSEA$Gene.Set.Name, levels = ResultsGSEA$Gene.Set.Name[order(as.numeric(ResultsGSEA$X..Genes.in.Gene.Set..K.), decreasing = FALSE)])
   ggplot(ResultsGSEA, aes_string(x = "X..Genes.in.Overlap..k.", y = "Gene.Set.Name", fill = "p.value")) +
     geom_col()+
     theme_bw() +
     #scale_fill_continuous(low="red", high="red4", name = "p.value", guide=guide_colorbar(reverse=TRUE)) +
     scale_fill_continuous(type = "gradient") +
     theme(axis.text = element_text(size = 12)) +
     # ggtitle(title) +
     xlab(NULL) + ylab(NULL)
 }
 
 
 
 ggplot(ResultsGSEA, aes_string(x = "X..Genes.in.Overlap..k.", y = "Gene.Set.Name", fill = "p.value")) +
   geom_col()+
   theme_bw() +
   scale_fill_continuous(low="red", high="red4", name = "p.value", guide=guide_colorbar(reverse=TRUE)) +
   #scale_fill_continuous(type = "gradient") +
   theme(axis.text = element_text(size = 12)) +
   # ggtitle(title) +
   xlab(NULL) + ylab(NULL)
 
 
 
 ggplot(ResultsGSEA, aes(x = "X..Genes.in.Overlap..k.", y = "Gene.Set.Name", fill = "p.value")) +
   geom_col()+
   theme_bw() +
   scale_fill_continuous(low="red", high="red4", name = "p.value", guide=guide_colorbar(reverse=TRUE)) +
   #scale_fill_continuous(type = "gradient") +
   theme(axis.text = element_text(size = 12)) +
   # ggtitle(title) +
   xlab(NULL) + ylab(NULL)
 
 
 
 ggplot(ResultsGSEA, aes(x = "X..Genes.in.Overlap..k.", y = "Gene.Set.Name")) +
   geom_bar(stat="identity", width=.5, fill="p.value") 
   
 
 
 
 
 
 #PCA tras el filtrado
 par(mfrow=c(1,2))
 plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
         col=pal[factor(targets$TISSUE)])
 legend("top", legend=levels(factor(targets$TISSUE)), text.col=pal,
        bg="white", cex=0.7)
 
 
 
 library(limma)
 
 plotMDS(mVals120, top=1000, gene.selection="common", 
         col=pal[factor(targets$TISSUE)])
 legend("top", legend=levels(factor(targets$TISSUE)), text.col=pal,
        bg="white", cex=0.7)
 
 
 tsne(mVals120,labels=as.factor(targets$TISSUE))

 
 
 
 plot(tsne_results$Y, col = "black", bg= pal[factor(targets$TISSUE)], pch = 21, cex = 1.5) 
 
 plot(tsne_results$Y, col = "black", bg= pal[factor(targets$TISSUE)], pch = 21, cex = 1.5)  +  geom_text_repel(label=targets$TISSUE) 
   
 
 tsneDF <- as.data.frame(tsne_results$Y)
 
 
 library(ggplot2)
 library(ggpubr)
 library(ggrepel)
 
  ggplot(tsneDF,aes(x = tsneDF$V1,y = tsneDF$V2)) +
   #labs(x=sonda) +
   geom_point() +
   geom_text_repel(label=targets$TISSUE) #+
   #geom_text_repel(label= ifelse(k$SAMPLE_NAME  %in% muestrasEnLab,k$SAMPLE_NAME,""), max.overlaps = Inf)  +
   #ggtitle(paste("Expression-methylation graph for the gene:",i,"that is found in:",variablePromotor))
  
  
  
  
  
  
  
  
  #filtramos lineas celulares
  CosmicSSEPIC <- CosmicCompleteCL[CosmicCompleteCL$SAMPLE_NAME %in% betas_resultEPIC$Sample.Name,]
  #Despues filtramos por genes. En mi caso obtengo 10360 lineas.
  CosmicSSEPICGenes <- CosmicSSEPIC[CosmicSSEPIC$GENE_NAME %in% joinGenes,]
  
  
  #CosmicSSEPICGenes <- CosmicCompleteCL[CosmicCompleteCL$GENE_NAME %in% listaGenesGencode,]
  
  
  
  pdf("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/graficas expresion-metilacion/GraficaEpicAllGencode2023.pdf")
  
  #Comentamos la lineas de UCSC o Gencode segun queramos usar cada una de las anotaciones. Ejecutamos 1 vez con cada anotacion y obtenemos graficas para todos los genes.
  #Tambien esta la posibilidad de ejecutar con todas las lineas celulares producto de mezclar 450K + EPIC o solo con EPIC.
  i <- "ACYP2"
  
  for (i in listaGenesGencode) {
    
    if(i %in% CosmicSSAll$GENE_NAME){ #Ejecuto el script solo si el gen esta en la lista de nombre de genes de cosmic
      
     # rm(sonda)
      
      #j <- assign(paste0("gen", i), CosmicSSEPICGenes[CosmicSSEPICGenes$GENE_NAME %in% i,])
      j <- assign(paste0("gen", i), CosmicSSAll[CosmicSSAll$GENE_NAME %in% i,])
      
      sonda <- DMPs0.15Signi[DMPs0.15Signi$GencodeBasicV12_NAME %in% i, "Name" ]
      
      sonda <- sonda[1]
      
      if(!is_empty(sonda) && sonda  %in% colnames(betasAll)){
        #expresionGen <- betas_resultEPIC[betas_resultEPIC$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
        expresionGen <- betasAll[betasAll$Sample.Name %in% j$SAMPLE_NAME ,c("Sample.Name", sonda)]
        
        #Si tengo mas de un valor para una linea celular hago la media.
        expresionGen <- as.data.table(expresionGen)
        expresionGen <- expresionGen[,lapply(.SD, mean), by = Sample.Name]
        
        k <- merge(j, expresionGen, by.x = "SAMPLE_NAME", by.y = "Sample.Name")
        
        #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$UCSC_RefGene_Name),"UCSC_RefGene_Group"]
        #variablePromotor <- DMPs0.15Signi[match(i,DMPs0.15Signi$GencodeBasicV12_NAME),"GencodeBasicV12_Group"]
        
        plot <- ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
          labs(x=sonda) +
          geom_point() +
          geom_smooth(method='lm') +
          stat_cor(method = "spearman" ) +
          theme_bw() +
          ggtitle(paste("Expression-methylation graph for the gene:",i,"that is found in:"))
        # ggtitle(paste("Expression-methylation graph for the gene:",i,"that is found in:",variablePromotor))
        
        print(plot)
        
      }
    }
  } 
  
  
  dev.off()
  
  print(plot)
  
  
  plot
  
  
  
  
 plot <-  ggplot(k,aes(x = k[ ,sonda],y = Z_SCORE)) +
    labs(x=sonda) +
    geom_point() +
    geom_smooth(method='lm') +
    stat_cor(method = "spearman" ) +
    theme_bw() +
     ggtitle(paste("Expression-methylation graph for the gene:",i,"that is found in:",variablePromotor))
  
   
plot 




for (i in listaGenesGencode) {
   CosmicSSAll[CosmicSSAll$GENE_NAME %in% i,"GENE_NAME"]
}





dim(DMPs0.15Signi[!DMPs0.15Signi$GencodeBasicV12_NAME %in% "",])
dim(DMPs0.15Signi[DMPs0.15Signi$GencodeBasicV12_NAME %in% "",])


length(BetaDiferenciaMedia0.15$Diferencia >= 0.15)

BetaDiferencia120 <- BetaDiferenciaMedia0.15[rownames(BetaDiferenciaMedia0.15) %in% rownames(DMPs0.15Signi),]

length(BetaDiferencia120$Diferencia <= -0.15)


dim(BetaDiferencia120[BetaDiferencia120$Diferencia <= 0.15,])

dim(BetaDiferencia120[BetaDiferencia120$Diferencia >= -0.15,])



DMPs0.15Signi66 <- DMPs0.15Signi[!DMPs0.15Signi$GencodeBasicV12_NAME %in% "",]

a <- DMPs0.15Signi66[DMPs0.15Signi66$Methyl450_Loci %in% TRUE,]


genDSTN$SAMPLE_NAME
write.table(genACYP2$SAMPLE_NAME, file="/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/listaLineasCelularesEpic+450k", sep=",", row.names=FALSE)


library(biomaRt)
mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",         
                            mart    = useMart("ENSEMBL_MART_ENSEMBL",       
                                              host    = "https://www.ensembl.org"))


geneSet <- c("ZNF608","IRF2BPL")
resultTable <- biomaRt::getBM(attributes = c("start_position","end_position","description"),       
                              filters    = "hgnc_symbol",       
                              values     = geneSet,         
                              mart       = mart)        

resultTable




DMPs0.15Signi <- readRDS("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/DMPs0.15Signi120.rds")


listaGenesGencode <- DMPs0.15Signi$GencodeBasicV12_NAME
listaGenesGencode <- na.omit(listaGenesGencode)
listaGenesGencode<- strsplit(as.character(listaGenesGencode), ";")
listaGenesGencode <-unlist(listaGenesGencode)

listaGenesGencode<-unique(listaGenesGencode)




lungFibrosis <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/dataAutopsias/GSE159655_Avg_Beta_Matrix.txt", TRUE, sep = "\t")


heatmap3(bValsfibrosis120, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 




heatmap3(bValsfibrosis105, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 


heatmap3(bValsfibrosis404, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 



heatmap3(bValscigar405, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 


library(EpiDISH)

data(centDHSbloodDMC.m)
#out.l <- epidish(bVals, centDHSbloodDMC.m, method = 'RPC')




#Utilizo funcion wRPC del paquete EPISCORE

targets <- read.metharray.sheet("/Users/jparra/Documents/ProyectoCovid/data/blueprint", pattern="SampleSheetAutop34.csv")

targetsAleix <- targets[!targets$EPIC.Sample.Name %in% c("L19-C-Lu","L22-C-Lu"),]

#bVals <- readRDS("/Users/jparra/Documents/ProyectoCovid/data/deconvolution/bValsAleix.rds")
Bvals_Patients <- read.csv("/Users/jparra/Documents/ProyectoCovid/data/deconvolution/Bvals_Patients.csv")

library(EpiSCORE)
rownames(Bvals_Patients)<- Bvals_Patients$ï..CpG_ID
Bvals_Patients <- Bvals_Patients[,2:33]

unique(rownames(Bvals_Patients))

#colnames(bVals) <- targets$EPIC.Sample.Name
#avDNAm.m <- constAvBetaTSS(Bvals_Patients,type="850k")

Bvals_PatientsM <- as.matrix(Bvals_Patients)
avDNAm.m <- constAvBetaTSS(Bvals_PatientsM,type="850k")

rownames(avDNAm.m)
unique(rownames(avDNAm.m))

#estF.o <- wRPC(avDNAm.m,ref=centDHSbloodDMC.m,useW=TRUE,wth=0.2,maxit=200);
estF.o <- wRPC(data=avDNAm.m,ref=centDHSbloodDMC.m,useW=F)

estF.o$estF
estF.o$ref

write.table(estF.o$estF, file="/Users/jparra/Documents/ProyectoCovid/data/blueprint/bigwig/celulasEpiProfile/deconvolution2.csv", sep=",", col.names=TRUE)




#epidish
estF.m <- epidish(Bvals_Patients,ref.m=centDHSbloodDMC.m,method="RPC",maxit=200)
#estF.m$ref

#estF.m$dataREF

estF.m$estF


colnames(estF.m$estF) <- c("Bcells","Natural Killer","CD4T","CD8T","Monocitos","Neutrofilos","Eosinofilos")





#pv <- wilcox.test(estF.o$estF[,"IC"] ~ phenoLUSC.lv$Cancer,alt="less")$p.value;
pv <- wilcox.test(estF.m$estF[,"B"] ~ targetsAleix$TISSUE,alt="less", exact=FALSE)
print(pv)



#La concentracion de celulas del sistema inmune es estadisticamente mas abundante en los pacientes con covid.

InmuneCellPercentage1<- as.data.frame(estF.o$estF[,"IC"])
InmuneCellPercentage1$group <- "control"
InmuneCellPercentage1$group[19:34] <- "COVID"
InmuneCellPercentage1$cellType <- "IC"
colnames(InmuneCellPercentage1) <- c("value", "group","cellType")

InmuneCellPercentage2<- as.data.frame(estF.o$estF[,"Fib"])
InmuneCellPercentage2$group <- "control"
InmuneCellPercentage2$group[19:34] <- "COVID"
InmuneCellPercentage2$cellType <- "Fib"
colnames(InmuneCellPercentage2) <- c("value", "group","cellType")


a <- rbind(InmuneCellPercentage1,InmuneCellPercentage2)

InmuneCellPercentage3<- as.data.frame(estF.o$estF[,"Endo"])
InmuneCellPercentage3$group <- "control"
InmuneCellPercentage3$group[19:34] <- "COVID"
InmuneCellPercentage3$cellType <- "Endo"
colnames(InmuneCellPercentage3) <- c("value", "group","cellType")

a <- rbind(a,InmuneCellPercentage3)


InmuneCellPercentage4<- as.data.frame(estF.o$estF[,"Epi"])
InmuneCellPercentage4$group <- "control"
InmuneCellPercentage4$group[19:34] <- "COVID"
InmuneCellPercentage4$cellType <- "Epi"
colnames(InmuneCellPercentage4) <- c("value", "group","cellType")

a <- rbind(a,InmuneCellPercentage4)


ggplot(a, aes(x = cellType, y = value, fill = group)) +    # Create boxplot chart in ggplot2
  geom_boxplot() + theme_bw()



ggplot(InmuneCellPercentage, aes(x = group, y = values, fill = group)) +    # Create boxplot chart in ggplot2
  geom_boxplot()




#Fisher Fibrosis BAL escalado o sin escalado
data <- c(7,7,0,30)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 8.956e-05

#Fisher E-cigar BAL con  escalado
data <- c(5,23,4,26)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.726

#Fisher E-cigar BAL sin escalado
data <- c(5,23,0,30)
data <- matrix(data, ncol=2)
fisher.test(data)
p-value = 0.02145

#epidish boold
estF.m <- epidish(Bvals_Patients,ref.m=centDHSbloodDMC.m,method="RPC",maxit=200)

estF.m$estF




ggplot(a, aes(x = cellType, y = value, fill = group)) +    # Create boxplot chart in ggplot2
  geom_boxplot() + theme_bw()


heatmap3(bValsfibrosis105, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 


library(gplots)
heatmap.2(bValsfibrosis105, ColSideColors=x, ColSideLabs="Signature",scale = "row",
         col= colorRampPalette(c("green", "black", "red"))(1024), method = "ward") 

library(ComplexHeatmap)
palette_greenred <- colorRampPalette( c("green", "black", "red"), space="rgb")(75)

Heatmap(bValsfibrosis405, name = "mat",cluster_rows=T,use_raster=TRUE,
        show_row_names=FALSE,show_column_names=FALSE,
        raster_device='CairoPNG', col = palette_greenred,
        row_dend_reorder = TRUE,show_row_dend = F, ColSideColors=y)


Heatmap(bValsfibrosis405, name = "mat",cluster_rows=T,use_raster=TRUE,
        show_row_names=FALSE,show_column_names=FALSE,
        raster_device='CairoPNG', col = palette_greenred,
        row_dend_reorder = TRUE,show_row_dend = F,bottom_annotation =anotacion,clustering_method_columns = "ward.D")

y

Heatmap(bValscigar405, name = "mat",cluster_rows=T,use_raster=TRUE,
        show_row_names=FALSE,show_column_names=FALSE,
        raster_device='CairoPNG', col = palette_greenred,
        row_dend_reorder = TRUE,show_row_dend = F,bottom_annotation =anotacion,clustering_method_columns = "ward.D")


clustering_method_columns = "ward.D",









#aparte


pheatmap(mat, annotation_col = anno)

pheatmap(mat, annotation_col = anno, fontsize_row = 6)

pheatmap(mat, annotation_col = anno, legend=T)

rownames(mat)

## Create the heatmap:
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(mat, annotation_col = anno)
setHook("grid.newpage", NULL, "replace")
#library(grid)
#grid.text("xlabel example", y=-0.07, gp=gpar(fontsize=16))
grid.text(rownames(mat), x=-0.07, rot=90, gp=gpar(fontsize=16))



pheatmap(mat, annotation_col = anno,cutree_rows=5)


pheatmap(mat, annotation_col = anno, fontsize_row = 12,cutree_rows=5, fontsize_col = 20)


par(mfrow = c(1, 2))
for (i in 1:length(colnames(fit.cont$contrasts))){
  plotMD(fit.cont, coef = i, status = summa.fit[, i], values = c(-1, 1))
  volcanoplot(fit.cont, coef = i, highlight = 100, 
              names = fit.cont$genes)
  
}


#16
res.list <- list()
for (i in 1:length(colnames(fit.cont$contrasts))){
  x <- topTable(fit.cont, coef = i, sort.by = "p", n = Inf, confint = T)
  res.list[[i]] <- x
  names(res.list)[i] <- colnames(fit.cont$contrasts)[i]
 # write.csv(x, paste0("results/", experiment.title, "_DEGs_", 
                      #colnames(fit.cont$contrasts)[i], ".csv"))
}

#19
library(qvalue)

#pacman::qvalue
synergy.pvalues <- res.list$Synergy$P.Value
pi1 <- 1 - qvalue(synergy.pvalues)$pi0
print(pi1)

plot.new()
text(0.4,0.75,labels=paste0("\n",round(pi1 * 100, 2), 
                            "% non-null \np-values and \n", 
                            round(sum(res.list$Synergy$adj.P.Val < 0.1) *
                                    100/length(res.list$Synergy$ensembl), 2), 
                            " % of genes with \nsynergy FDR < 0.1"))


print(h)
plotSA(fit.cont, main = "Final model: Mean-variance trend", ylab = "Sqrt( standard deviation )")


vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

  


res.list[[5]] %>%
  ggplot(aes(x = logFC,
             y = -log10(P.Value),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  geom_label_repel(data = sig_up_il_genes, # Add labels last to appear as the top layer  
                   aes(label = ID),
                   force = 2,
                   nudge_y = 1) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 5)),       
                     limits = c(-10, 10)) 




res.list[[5]] %>%
  ggplot(aes(x = logFC,
             y = -log10(P.Value),
             #fill = gene_type,    
             #size = gene_type,
             #alpha = gene_type
             label= 
             )) + 
  geom_point(aes(colour = gene_type), 
             alpha = 0.2, 
             shape = 16,
             size = 1) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  geom_label_repel(data = sig_up_il_genes, # Add labels last to appear as the top layer  
                   aes(label = ID,colour = gene_type),
                   force = 2,
                   nudge_y = 1) +
  #scale_fill_manual(values = cols) + # Modify point colour
  #scale_size_manual(values = sizes) + # Modify point size
  #scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 5)),       
                     limits = c(-10, 10))


res.list[[5]] %>%
  ggplot(aes(x = logFC,
             y = -log10(P.Value),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type
             #,label=ID
  )) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  geom_label_repel(data = sig_up_il_genes, # Add labels last to appear as the top layer  
                   aes(label = ID),
                  # force = 2,
                  # nudge_y = 1,
                   max.overlaps = Inf) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 5)),       
                     limits = c(-10, 10)) 




res.list[[5]] %>%
  ggplot(aes(x = logFC,
             y = -log10(P.Value),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type
             # ,label=ID
  )) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-5, 5),
             linetype = "dashed") +
  geom_label_repel(data = sig_up_il_genes, # Add labels last to appear as the top layer  
                   aes(label = ID,colour = gene_type),
                   force = 2,
                   nudge_y = 1,max.overlaps = Inf) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 5)),       
                     limits = c(-20, 20)) 




pheatmap(mat, annotation_col = anno)
plot <-pheatmap(mat, annotation_col = anno, fontsize_row = 12,cutree_rows=12, fontsize_col = 20)
listAttributes(mart)

mart@attributes$name[mart@attributes$name %in% "ensembl_gene_i"]
plot
pheatmap(mat[,1:27], annotation_col = anno)

pheatmap(heatmapCounts, annotation_col = anno)
pheatmap(mat, annotation_col = anno)
pheatmap(heatmapCountsTOPUnique[,1:12], annotation_col = anno)




pheatmap(tmp,
         #kmeans_k = 30, 
         #annotation_row = rownames(tmp),
         #cellwidth = 70, cellheight = 5,
         border_color = NA, 
         breaks=breaks,
         cluster_cols = F,
         show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(117),
         main = paste0("logFC expected vs. measured:\n",
                       levels(log2FC.matrixOrderUnique$magnitude.syn)[i]))
print(h)



for (i in 2:5){
  breaks <- c(seq(-6, -0.3,by=0.1),seq(0.3, 6,by=0.1))
  breaks <- append(breaks, -9,0)
  breaks <- append(breaks, 9)
  tmp <- log2FC.matrixOrderUniqueG2[log2FC.matrixOrderUniqueG2$magnitude.syn == 
                                      levels(log2FC.matrixOrderUniqueG2$magnitude.syn)[i],
                                    c("Additive.logFC","Combinatorial.logFC")]
  h <- pheatmap(tmp,
                #kmeans_k = 30, 
                #annotation_row = rownames(tmp),
                #cellwidth = 70, cellheight = 5,
                border_color = NA, 
                breaks=breaks,
                cluster_cols = F,
                show_rownames = T,
                color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(117),
                main = paste0("logFC expected vs. measured:\n",
                              levels(log2FC.matrixOrderUniqueG2$magnitude.syn)[i]))
  print(h)
}


for (i in 1:5){
  breaks <- c(seq(-6, -0.3,by=0.1),seq(0.3, 6,by=0.1))
  breaks <- append(breaks, -9,0)
  breaks <- append(breaks, 9)
  tmp <- log2FC.matrixg2cellCycle[log2FC.matrixg2cellCycle$magnitude.syn == 
                                    levels(log2FC.matrixg2cellCycle$magnitude.syn)[i],
                                  c("Additive.logFC","Combinatorial.logFC")]
  if(!dim(tmp)[[1]]==0 )
    h <- pheatmap(tmp,
                  #kmeans_k = 30, 
                  #annotation_row = rownames(tmp),
                  #cellwidth = 70, cellheight = 5,
                  fontsize = 6,
                  border_color = NA, 
                  #breaks=breaks,
                  cluster_cols = F,
                  cluster_rows = F,
                  show_rownames = T,
                  color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(117),
                  main = paste0("logFC expected vs. measured:\n",
                                levels(log2FC.matrixg2cellCycle$magnitude.syn)[i]))
  print(h)
}



pheatmap(mat, annotation_col = anno)




listAttributes(mart)


install.packages("tidyverse")

install.packages("dplyr")

install.packages("Rtools")
