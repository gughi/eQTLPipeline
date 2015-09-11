#######################
### Sentinalisation ###
#######################  
rm(list=ls())

library(devtools)
library(doParallel)
library(foreach)
load_all()

setwd("/home/guelfi/eQTLPipeline/")
cl <- makeCluster(15)
clusterExport(cl,c("LDsentinalisationIntergenic","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()
## path we put the results
pathFinalSentinalised <-"data/results/intergenic/resMatrixEQTL/sentinalised/PUTM"
## path where get the unsentinalised eQTLs
pathUnsentinalised <- "data/results/intergenic/resMatrixEQTL/PUTM/"



regions <- read.delim(pipe(paste0("ls ",pathUnsentinalised)),header=F)

load("data/expr/normalisedCounts/intergenic/resids.PUTM.rda")
load("data/general/sampleInfo.rda")

IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

pathFinalSentinalised <-"data/results/intergenic/resMatrixEQTL/sentinalised/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)
pathFinalSentinalised <-"data/results/intergenic/resMatrixEQTL/sentinalised/PUTM/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)

dir.create(file.path("tmp", "PUTM"),showWarnings=FALSE)

# load the genetic PCs
my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)


Sys.time()
foreach(i=1:nrow(regions))%dopar%LDsentinalisationIntergenic(exprIntergenic=resids,
                                                   regID=regions[i,1],
                                                        pathFinalSentinalised=pathFinalSentinalised,
                                                        pathUnsentinalised=pathUnsentinalised,
                                                        FDRthr=0.10,
                                                        my.covTMP=my.covTMP,
                                                        snpLocation="data/snps/byRegion/PUTM/",
                                                        tmpFolder="tmp/")
Sys.time()
stopCluster(cl)

rm(list=ls())


cl <- makeCluster(15)
clusterExport(cl,c("LDsentinalisationIntergenic","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()
## path we put the results
pathFinalSentinalised <-"data/results/intergenic/resMatrixEQTL/sentinalised/SNIG"
## path where get the unsentinalised eQTLs
pathUnsentinalised <- "data/results/intergenic/resMatrixEQTL/SNIG/"



regions <- read.delim(pipe(paste0("ls ",pathUnsentinalised)),header=F)

load("data/expr/normalisedCounts/intergenic/resids.SNIG.rda")
load("data/general/sampleInfo.rda")

IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

pathFinalSentinalised <-"data/results/intergenic/resMatrixEQTL/sentinalised/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)
pathFinalSentinalised <-"data/results/intergenic/resMatrixEQTL/sentinalised/SNIG/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)

dir.create(file.path("tmp", "SNIG"),showWarnings=FALSE)


# load the genetic PCs
my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)


Sys.time()
foreach(i=1:nrow(regions))%dopar%LDsentinalisationIntergenic(exprIntergenic=resids,
                                                             regID=regions[i,1],
                                                             pathFinalSentinalised=pathFinalSentinalised,
                                                             pathUnsentinalised=pathUnsentinalised,
                                                             FDRthr=0.10,
                                                             my.covTMP=my.covTMP,
                                                             snpLocation="data/snps/byRegion/SNIG/",
                                                             tmpFolder="tmp/")
Sys.time()
stopCluster(cl)











