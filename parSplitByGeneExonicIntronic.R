#########################################################################
### We split the expression by gene using the parallelisation of the jobs
#########################################################################

# delete all the objects
rm(list=ls())

library(devtools)
library(doParallel)
library(foreach)
load_all()
## detectCores()
## [1] 24

setwd("/home/guelfi/eQTLPipeline/")
cl <- makeCluster(15)
clusterExport(cl,"splitExprByGene")
registerDoParallel(cl)
getDoParWorkers()
## [1] 16
## we load the expression data this to save the loading time
ensemblGenes <- read.delim(file="data/general/ensemblGenes.txt", as.is=T,header=T)
ensemblRef <- read.delim(file="data/general/ensemblRef.txt", as.is=T,header=T)
pathResExpr <- "data/expr/normalisedCounts/genic/ExonIntrons/"
## load the data from PUTM
load("data/expr/normalisedCounts/genic/ExonIntrons/resids.PUTM.rda")
exprPUTM <- resids
## load the data from SNIG
rm(resids)
load("data/expr/normalisedCounts/genic/ExonIntrons/resids.SNIG.rda")
exprSNIG <- resids
rm(resids)
## load the covariates
load("data/expr/normalisedCounts/genic/ExonIntrons/RPKM.cqn.PUTM")
rm(covs,RPKM.cqn)

load("data/expr/normalisedCounts/genic/ExonIntrons/RPKM.cqn.SNIG")
rm(covs,RPKM.cqn)

Sys.time()
foreach(i=1:nrow(ensemblGenes))%dopar%splitExprByGene(i,ensemblRef,ensemblGenes,PUTM,exprPUTM,SNIG,exprSNIG,pathResExpr)
Sys.time()
stopCluster(cl)
rm(cl)

