############################
## Run the eQTL analysis ###
############################

# delete all the objects

library(devtools)
library(doParallel)
library(foreach)
load_all()


# Caprica
setwd("/home/seb/projectsR/eQTLPipeline/")
## load the annotation for the exon-exon junctions
load("data/expr/rawCounts/genic/fullExExJun.rda")

# Apollo
# setwd("/home/guelfi/eQTLPipeline/")
cl <- makeCluster(20)
clusterExport(cl,c("runCisEQTLExExJun","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()

## load the residual corrected expression
load("data/expr/normalisedCounts/genic/exonExonJunc/resids.PUTM.rda")
## load the information about tissue
load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.PUTM")
## remove the object we don't need
rm(covs,RPKM.cqn)


load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(IDs)

## Apollo
## snpLocation <- "/home/guelfi/eQTL/snps/byGene/"
## Caprica
snpLocation <- "/home/seb/eQTL/snps/byGene/"

# Apollo
# outputFolder <- "/home/guelfi/eQTLPipeline/data/results/genic/exonExonJunc/resMatrixEQTL/PUTM/"
# fullResults <- "/home/guelfi/eQTLPipeline/data/results/genic/exonExonJunc/fullResults/PUTM/"

# Caprica
outputFolder <- "/home/seb/projectsR/eQTLPipeline/data/results/genic/exonExonJunc/resMatrixEQTL/PUTM/"
fullResults <- "/home/seb/projectsR/eQTLPipeline/data/results/genic/exonExonJunc/fullResults/PUTM/"

# Apollo
# my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)
# Caprica
my.covTMP <- read.table.rows(paste0("/home/seb/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)

Sys.time()
foreach(i=1:ncol(resids))%dopar%runCisEQTLExExJun(i=i,resids=resids,
                                                     mapExon=mapExon,
                                                     snpLocation=snpLocation,
                                                     outputFolder=outputFolder,
                                                     my.covTMP=my.covTMP,
                                                     fullResults=fullResults)
Sys.time()
stopCluster(cl)
rm(cl)


library(devtools)
library(doParallel)
library(foreach)
load_all()


# Caprica
setwd("/home/seb/projectsR/eQTLPipeline/")
## load the annotation for the exon-exon junctions
load("data/expr/rawCounts/genic/fullExExJun.rda")

# Apollo
# setwd("/home/guelfi/eQTLPipeline/")
cl <- makeCluster(20)
clusterExport(cl,c("runCisEQTLExExJun","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()



## load the residual corrected expression
load("data/expr/normalisedCounts/genic/exonExonJunc/resids.SNIG.rda")
## load the information about tissue
load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.SNIG")
## remove the object we don't need
rm(covs,RPKM.cqn)

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(IDs)

## Apollo
## snpLocation <- "/home/guelfi/eQTL/snps/byGene/"
## Caprica
snpLocation <- "/home/seb/eQTL/snps/byGene/"

# Apollo
# outputFolder <- "/home/guelfi/eQTLPipeline/data/results/genic/exonExonJunc/resMatrixEQTL/SNIG/"
# fullResults <- "/home/guelfi/eQTLPipeline/data/results/genic/exonExonJunc/fullResults/SNIG/"

# Caprica
outputFolder <- "/home/seb/projectsR/eQTLPipeline/data/results/genic/exonExonJunc/resMatrixEQTL/SNIG/"
fullResults <- "/home/seb/projectsR/eQTLPipeline/data/results/genic/exonExonJunc/fullResults/SNIG/"

# Apollo
# my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)
# Caprica
my.covTMP <- read.table.rows(paste0("/home/seb/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)

Sys.time()
foreach(i=1:ncol(resids))%dopar%runCisEQTLExExJun(i=i,resids=resids,
                                                  mapExon=mapExon,
                                                  snpLocation=snpLocation,
                                                  outputFolder=outputFolder,
                                                  my.covTMP=my.covTMP,
                                                  fullResults=fullResults)
Sys.time()
stopCluster(cl)
rm(cl)



