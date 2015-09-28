############################
## Run the eQTL analysis ###
############################

# delete all the objects
    
  library(devtools)
  library(doParallel)
  library(foreach)
  load_all()

  setwd("/home/guelfi/eQTLPipeline/")
  cl <- makeCluster(7)
  clusterExport(cl,c("runCisEQTLGene","read.table.rows"))
  registerDoParallel(cl)
  getDoParWorkers()
 
  ## load the residual corrected expression
  load("data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda")
  ## load the information about tissue
  load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
  ## remove the object we don't need
  rm(covs,RPKM.cqn)


  load("data/general/sampleInfo.rda")
  IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
  IDs <- gsub("/","_",IDs)
  rownames(resids) <- IDs
  rm(IDs)
 
 snpLocation <- "/home/guelfi/eQTL/snps/byGene/"
 
 outputFolder <- "/home/guelfi/eQTLPipeline/data/results/genic/geneIntronic/resMatrixEQTL/"
 dir.create(file.path(outputFolder),showWarnings=FALSE)
 outputFolder <- "/home/guelfi/eQTLPipeline/data/results/genic/geneIntronic/resMatrixEQTL/PUTM/"
 
 fullResults <- "/home/guelfi/eQTLPipeline/data/results/genic/geneIntronic/fullResults/"
 dir.create(file.path(fullResults),showWarnings=FALSE)
 fullResults <- "/home/guelfi/eQTLPipeline/data/results/genic/geneIntronic/fullResults/PUTM/"

 
 
 
 my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)
 
 Sys.time()
 foreach(i=1:ncol(resids))%dopar%runCisEQTLGene(i=i,resids=resids,
                                                     snpLocation=snpLocation,
                                                     outputFolder=outputFolder,
                                                     my.covTMP=my.covTMP,
                                                     fullResults=fullResults)

Sys.time()
stopCluster(cl)
rm(cl)


# rm(list=ls())
# 
# cl <- makeCluster(7)
# clusterExport(cl,c("runCisEQTLExons","read.table.rows"))
# registerDoParallel(cl)
# getDoParWorkers()
# 
# 
# ## load the residual corrected expression
# load("data/expr/normalisedCounts/genic/exons/resids.SNIG.rda")
# ## load the information about tissue
# load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.SNIG")
# ## remove the object we don't need
# rm(covs,RPKM.cqn)
# 
# 
# load("data/general/sampleInfo.rda")
# IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
# IDs <- gsub("/","_",IDs)
# rownames(resids) <- IDs
# rm(IDs)
# 
# snpLocation <- "/home/guelfi/eQTL/snps/byGene/"
# 
# outputFolder <- "/home/guelfi/eQTLPipeline/data/results/genic/exons/resMatrixEQTL/SNIG/"
# fullResults <- "/home/guelfi/eQTLPipeline/data/results/genic/exons/fullResults/SNIG/"
# 
# my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)
# 
# Sys.time()
# foreach(i=1:ncol(resids))%dopar%runCisEQTLExons(i=i,resids=resids,
#                                                 snpLocation=snpLocation,
#                                                 outputFolder=outputFolder,
#                                                 my.covTMP=my.covTMP,
#                                                 fullResults=fullResults)
# 
# Sys.time()
# stopCluster(cl)
# rm(cl)



