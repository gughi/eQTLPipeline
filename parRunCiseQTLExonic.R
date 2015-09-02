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
cl <- makeCluster(7)
clusterExport(cl,"runCisEQTL","read.table.rows")
registerDoParallel(cl)
getDoParWorkers()


ensemblGenes <- read.delim(pipe("ls data/expr/normalisedCounts/genic/geneExons/byGene_snps1Mb/"),header=F)
numGenes <- nrow(ensemblGenes)

Sys.time()
foreach(i=1:numGenes)%dopar%runCisEQTL(i=i,ensemblGenes=ensemblGenes,
                                       exprLocation="expr/normalisedCounts/genic/geneExons/byGene_snps1Mb/",
                                       snpLocation="/home/seb/eQTL/snps/byGene/",
                                       outputFolder="data/results/genic/geneExons/resMatrixEQTL",
                                       genotypeFile="/home/seb/plinkOutput/eigenvec",
                                       fullResults="data/results/genic/geneExons/fullResults/")
##foreach(i=1:30)%dopar%runCisEQTL(i,ensemblGenes)
Sys.time()
stopCluster(cl)
rm(cl)
