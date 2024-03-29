#######################
### Sentinalisation ###
#######################  
rm(list=ls())

library(devtools)
library(doParallel)
library(foreach)
load_all()

setwd("/home/guelfi/eQTLPipeline/")
cl <- makeCluster(7)
clusterExport(cl,c("LDsentinalisation2","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()
## path where get the unsentinalised eQTLs
pathUnsentinalised <- "data/results/genic/geneIntronic/resMatrixEQTL/PUTM/"

genes <- read.delim(pipe(paste0("ls ",pathUnsentinalised)),header=F)


## load the residual corrected expression
load("data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda")

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

pathFinalSentinalised <-"data/results/genic/geneIntronic/resMatrixEQTL/sentinalised/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)
pathFinalSentinalised <-"data/results/genic/geneIntronic/resMatrixEQTL/sentinalised/PUTM/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)

dir.create(file.path("tmp", "PUTM"),showWarnings=FALSE)

# load the genetic PCs
my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)


Sys.time()
foreach(i=1:nrow(genes))%dopar%LDsentinalisation2(resids=resids,
                                                             regID=genes[i,1],
                                                             pathFinalSentinalised=pathFinalSentinalised,
                                                             pathUnsentinalised=pathUnsentinalised,
                                                             FDRthr=0.10,
                                                             my.covTMP=my.covTMP,
                                                             snpLocation="/home/guelfi/eQTL/snps/byGene/",
                                                             tmpFolder="tmp/")
Sys.time()
stopCluster(cl)

rm(list=ls())

cl <- makeCluster(7)
clusterExport(cl,c("LDsentinalisation2","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()
## path where get the unsentinalised eQTLs
pathUnsentinalised <- "data/results/genic/geneIntronic/resMatrixEQTL/SNIG/"

genes <- read.delim(pipe(paste0("ls ",pathUnsentinalised)),header=F)


## load the residual corrected expression
load("data/expr/normalisedCounts/genic/geneIntronic/resids.SNIG.rda")

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

pathFinalSentinalised <-"data/results/genic/geneIntronic/resMatrixEQTL/sentinalised/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)
pathFinalSentinalised <-"data/results/genic/geneIntronic/resMatrixEQTL/sentinalised/SNIG/"
dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)

dir.create(file.path("tmp", "SNIG"),showWarnings=FALSE)

# load the genetic PCs
my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)


Sys.time()
foreach(i=1:nrow(genes))%dopar%LDsentinalisation2(resids=resids,
                                                    regID=genes[i,1],
                                                    pathFinalSentinalised=pathFinalSentinalised,
                                                    pathUnsentinalised=pathUnsentinalised,
                                                    FDRthr=0.10,
                                                    my.covTMP=my.covTMP,
                                                    snpLocation="/home/guelfi/eQTL/snps/byGene/",
                                                    tmpFolder="tmp/")
Sys.time()
stopCluster(cl)


