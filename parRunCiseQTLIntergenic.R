############################
## Run the eQTL analysis ###
############################

# delete all the objects
rm(list=ls())

library(devtools)
library(doParallel)
library(foreach)
load_all()

setwd("/home/guelfi/eQTLPipeline/")
cl <- makeCluster(15)
clusterExport(cl,c("runCisEQTLIntergenic","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()



load("data/general/sampleInfo.rda")
load("data/expr/normalisedCounts/intergenic/resids.PUTM.rda")

IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

snpLocation <- "data/snps/byRegion/PUTM/"
outputFolder <- "data/results/intergenic/"
dir.create(file.path(outputFolder),showWarnings=FALSE)
dir.create(file.path(outputFolder,"resMatrixEQTL/"),showWarnings=FALSE)
dir.create(file.path(paste0(outputFolder,"resMatrixEQTL/"), "PUTM"),showWarnings=FALSE)
outputFolder <- paste0("data/results/intergenic/resMatrixEQTL/PUTM/")

fullResults <- "data/results/intergenic/fullResults/"
dir.create(file.path(fullResults),showWarnings=FALSE)
dir.create(file.path(paste0(fullResults), "PUTM"),showWarnings=FALSE)
fullResults <- "data/results/intergenic/fullResults/PUTM/"

my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)


Sys.time()
foreach(i=1:ncol(resids))%dopar%runCisEQTLIntergenic(i=i,exprIntergenic=resids,
                                       snpLocation="data/snps/byRegion/PUTM/",
                                       outputFolder=outputFolder,
                                       my.covTMP=my.covTMP,
                                       fullResults=fullResults)

Sys.time()
stopCluster(cl)
rm(cl)


# delete all the objects
rm(list=ls())

cl <- makeCluster(15)
clusterExport(cl,c("runCisEQTLIntergenic","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()


load("data/general/sampleInfo.rda")
load("data/expr/normalisedCounts/intergenic/resids.SNIG.rda")

IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

snpLocation <- "data/snps/byRegion/SNIG/"
outputFolder <- "data/results/intergenic/"
dir.create(file.path(outputFolder),showWarnings=FALSE)
dir.create(file.path(outputFolder,"resMatrixEQTL/"),showWarnings=FALSE)
dir.create(file.path(paste0(outputFolder,"resMatrixEQTL/"), "SNIG"),showWarnings=FALSE)
outputFolder <- paste0("data/results/intergenic/resMatrixEQTL/SNIG/")

fullResults <- "data/results/intergenic/fullResults/"
dir.create(file.path(fullResults),showWarnings=FALSE)
dir.create(file.path(paste0(fullResults), "SNIG"),showWarnings=FALSE)
fullResults <- "data/results/intergenic/fullResults/SNIG/"

my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)


Sys.time()
foreach(i=1:ncol(resids))%dopar%runCisEQTLIntergenic(i=i,exprIntergenic=resids,
                                                     snpLocation="data/snps/byRegion/SNIG/",
                                                     outputFolder=outputFolder,
                                                     my.covTMP=my.covTMP,
                                                     fullResults=fullResults)

Sys.time()
stopCluster(cl)
rm(cl)




