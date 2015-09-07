
setwd("/home/guelfi/eQTLPipeline")
sink("logExonicIntornic.log")
nCores <- 15
cat(paste("Number of cores",nCores,"\n"))
library(devtools)
load_all()

cat(paste("Processing PUTM region \n"))

## Intergenic data was counts were generated using DerFinder
load("data/expr/rawCounts/intergenic/PUTM/allChromosomes.rda")
## load data has 2 different objects: annotation,coverage
## load the sample info
load("data/general/sampleInfo.rda")
## we select only samples from PUTM

PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
rm(sampleInfo)

exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
colnames(exprIntergenic) <- c("chr","start","end","width",colnames(coverage))
rm(annotation,coverage)    


## select regions with length => 100bp
exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
# We don't do filtering since it doesn't make any since for DERFINDER; 
# regions are detected based on the expression
# creation of identifiers for the regions
IDs <- paste0("DER",c(1:nrow(exprIntergenic)))
rownames(exprIntergenic) <- IDs
rm(IDs)

cat(paste("number of intergenic regions >100bp",nrow(exprIntergenic)))



# load the library size
librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content intergenic and region length        
length <- exprIntergenic$width    
names(length) <- rownames(exprIntergenic)
stopifnot(identical(colnames(as.matrix(exprIntergenic[,as.character(PUTM$A.CEL_file)]))
                    ,names(librarySize)))
stopifnot(identical(rownames(exprIntergenic),names(length)))              

#convert in RPKM
RPKM.std <- RPKM(as.matrix(exprIntergenic[,as.character(PUTM$A.CEL_file)])
                 , NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std) 
rm(RPKM.std,length)
exprIntergenic <- exprIntergenic[as.character(genesList),]

cat(paste("number of regions included in the analysis:",nrow(exprIntergenic)))

cat("Calculating the GC content")
GCcontent <- GCcalculation(exprIntergenic[,1:4],genRef="/home/ukbec/bowtie2Index/genome37.72.fa"
                           ,pathBedtools = "/apps/BEDTools/2.24.0/bin/bedtools")
GCcontent <- as.data.frame(cbind(GCcontent[as.character(rownames(exprIntergenic)),],
                                 exprIntergenic$width))    
rownames(GCcontent) <- rownames(exprIntergenic)
colnames(GCcontent) <- c("GCcontent","length")



library(cqn)
library(scales)

stopifnot(identical(colnames(exprIntergenic[,as.character(PUTM$A.CEL_file)]),names(librarySize)))
stopifnot(identical(rownames(exprIntergenic[,as.character(PUTM$A.CEL_file)]),rownames(GCcontent)))  

cat("Conditional quantile normalisation")

my.cqn <- cqn(exprIntergenic[,as.character(PUTM$A.CEL_file)], lengths = GCcontent$length,x = GCcontent$GCcontent, sizeFactors = librarySize , verbose = TRUE)


png(paste0("plots/intergenic/CQNPUTM.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset


PUTM$U.Region_simplified <- NULL
covs <- PUTM 
rownames(covs) <- covs$A.CEL_file
#convert the female and male info in numeric
covs[covs=="M"]=0
covs[covs=="F"]=1
covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
covs[,5] <- as.numeric(covs[,5])
covs[,6] <- as.numeric(covs[,6])
lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
rownames(lanes) <- gsub("CEL","",rownames(lanes))
covs <- cbind(covs,librarySize[as.character(rownames(covs))])
covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                    "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")


starStopReg <- exprIntergenic[,1:4]

save(RPKM.cqn,PUTM,covs,starStopReg,file="data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")

  rm(list=ls())
  
  cat("Residual correction")
  
  load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
  
  ##doSwamp(RPKM.cqn,covs)
  
  PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
  PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]
  
  resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                                 "data/expr/normalisedCounts/intergenic/resids.PUTM.rda")
  





############
### SNIG ###
############

nCores<-15
load_all()

cat(paste("Processing SNIG region \n"))

## Intergenic data was counts were generated using DerFinder
load("data/expr/rawCounts/intergenic/SNIG/allChromosomes.rda")
## load data has 2 different objects: annotation,coverage
## load the sample info
load("data/general/sampleInfo.rda")
## we select only samples from SNIG

SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),]
rm(sampleInfo)

exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
colnames(exprIntergenic) <- c("chr","start","end","width",colnames(coverage))
rm(annotation,coverage)    


## select regions with length => 100bp
exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
# We don't do filtering since it doesn't make any since for DERFINDER; 
# regions are detected based on the expression
# creation of identifiers for the regions
IDs <- paste0("DER",c(1:nrow(exprIntergenic)))
rownames(exprIntergenic) <- IDs
rm(IDs)

cat(paste("number of intergenic regions >100bp",nrow(exprIntergenic)))



# load the library size
librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(SNIG$A.CEL_file),]
names(librarySize) <- as.character(SNIG$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content intergenic and region length        
length <- exprIntergenic$width    
names(length) <- rownames(exprIntergenic)
stopifnot(identical(colnames(as.matrix(exprIntergenic[,as.character(SNIG$A.CEL_file)]))
                    ,names(librarySize)))
stopifnot(identical(rownames(exprIntergenic),names(length)))              

#convert in RPKM
RPKM.std <- RPKM(as.matrix(exprIntergenic[,as.character(SNIG$A.CEL_file)])
                 , NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std) 
rm(RPKM.std,length)
exprIntergenic <- exprIntergenic[as.character(genesList),]

cat(paste("number of regions included in the analysis:",nrow(exprIntergenic)))

cat("Calculating the GC content")
GCcontent <- GCcalculation(exprIntergenic[,1:4],genRef="/home/ukbec/bowtie2Index/genome37.72.fa"
                           ,pathBedtools = "/apps/BEDTools/2.24.0/bin/bedtools")
GCcontent <- as.data.frame(cbind(GCcontent[as.character(rownames(exprIntergenic)),],
                                 exprIntergenic$width))    
rownames(GCcontent) <- rownames(exprIntergenic)
colnames(GCcontent) <- c("GCcontent","length")



library(cqn)
library(scales)

stopifnot(identical(colnames(exprIntergenic[,as.character(SNIG$A.CEL_file)]),names(librarySize)))
stopifnot(identical(rownames(exprIntergenic[,as.character(SNIG$A.CEL_file)]),rownames(GCcontent)))  

cat("Conditional quantile normalisation")

my.cqn <- cqn(exprIntergenic[,as.character(SNIG$A.CEL_file)], lengths = GCcontent$length,x = GCcontent$GCcontent, sizeFactors = librarySize , verbose = TRUE)


png(paste0("plots/intergenic/CQNSNIG.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset


SNIG$U.Region_simplified <- NULL
covs <- SNIG 
rownames(covs) <- covs$A.CEL_file
#convert the female and male info in numeric
covs[covs=="M"]=0
covs[covs=="F"]=1
covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
covs[,5] <- as.numeric(covs[,5])
covs[,6] <- as.numeric(covs[,6])
lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
rownames(lanes) <- gsub("CEL","",rownames(lanes))
covs <- cbind(covs,librarySize[as.character(rownames(covs))])
covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                    "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")


save(RPKM.cqn,SNIG,covs,starStopReg,file="data/expr/normalisedCounts/intergenic/RPKM.cqn.SNIG")


cat("Residual correction")

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.SNIG")

##doSwamp(RPKM.cqn,covs)

PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                               "data/expr/normalisedCounts/intergenic/resids.SNIG.rda")


### split by region

rm(list=ls())

nCores <- 15

setwd("/home/guelfi/eQTLPipeline")
cat("starting splitting the expression \n")
writeSH(nameSH="splitByRegion.sh",logName="splitByRegion",
        cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/parSplitByGeneExonicIntronic.R"),numThreads=(nCores+1))




