#######################
### Sentinalisation ###
#######################  
rm(list=ls())

library(devtools)
library(doParallel)
library(foreach)
load_all()

# Caprica
setwd("/home/seb/projectsR/eQTLPipeline/")
# Apollo
# setwd("/home/guelfi/eQTLPipeline/")

cl <- makeCluster(20)
clusterExport(cl,c("LDsentExonExonJunc","read.table.rows"))
registerDoParallel(cl)

# getDoParWorkers()
# ## path where get the unsentinalised eQTLs
pathUnsentinalised <- "data/results/genic/exonExonJunc/resMatrixEQTL/PUTM/"
junctions <- read.delim(pipe(paste0("ls ",pathUnsentinalised)),header=F)

## load the residual corrected expression
load("data/expr/normalisedCounts/genic/exonExonJunc/resids.PUTM.rda")
 
load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)
 
pathFinalSentinalised <-"data/results/genic/exonExonJunc/resMatrixEQTL/sentinalised/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)
pathFinalSentinalised <-"data/results/genic/exonExonJunc/resMatrixEQTL/sentinalised/PUTM/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)

dir.create(file.path("tmp", "PUTM"),showWarnings=FALSE)
 
# load the genetic PCs
my.covTMP <- read.table.rows(paste0("/home/seb/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)

## load mapping file
load("data/expr/rawCounts/genic/fullExExJun.rda")
rm(map,expr)

Sys.time()
foreach(i=1:nrow(junctions))%dopar%LDsentExonExonJunc(resids=resids,
                                                    regID=junctions[i,1],
                                                    mapExon=mapExon,
                                                    pathFinalSentinalised=pathFinalSentinalised,
                                                    pathUnsentinalised=pathUnsentinalised,
                                                    FDRthr=0.10,
                                                    my.covTMP=my.covTMP,
                                                    snpLocation="/home/seb/eQTL/snps/byGene/",
                                                    tmpFolder="tmp/")
Sys.time()
stopCluster(cl)

rm(list=ls())


# Caprica
setwd("/home/seb/projectsR/eQTLPipeline/")
# Apollo
# setwd("/home/guelfi/eQTLPipeline/")

cl <- makeCluster(20)
clusterExport(cl,c("LDsentExonExonJunc","read.table.rows"))
registerDoParallel(cl)

# getDoParWorkers()
# ## path where get the unsentinalised eQTLs
pathUnsentinalised <- "data/results/genic/exonExonJunc/resMatrixEQTL/SNIG/"
junctions <- read.delim(pipe(paste0("ls ",pathUnsentinalised)),header=F)

## load the residual corrected expression
load("data/expr/normalisedCounts/genic/exonExonJunc/resids.SNIG.rda")

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

pathFinalSentinalised <-"data/results/genic/exonExonJunc/resMatrixEQTL/sentinalised/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)
pathFinalSentinalised <-"data/results/genic/exonExonJunc/resMatrixEQTL/sentinalised/SNIG/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)

dir.create(file.path("tmp", "SNIG"),showWarnings=FALSE)

# load the genetic PCs
my.covTMP <- read.table.rows(paste0("/home/seb/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)

## load mapping file
load("data/expr/rawCounts/genic/fullExExJun.rda")
rm(map,expr)

Sys.time()
foreach(i=1:nrow(junctions))%dopar%LDsentExonExonJunc(resids=resids,
                                                      regID=junctions[i,1],
                                                      mapExon=mapExon,
                                                      pathFinalSentinalised=pathFinalSentinalised,
                                                      pathUnsentinalised=pathUnsentinalised,
                                                      FDRthr=0.10,
                                                      my.covTMP=my.covTMP,
                                                      snpLocation="/home/seb/eQTL/snps/byGene/",
                                                      tmpFolder="tmp/")
Sys.time()
stopCluster(cl)


