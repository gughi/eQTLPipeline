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
clusterExport(cl,c("runCisEQTL","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()


ensemblGenes <- read.delim(pipe("ls data/expr/normalisedCounts/genic/geneExons/byGene_snps1Mb/"),header=F)
numGenes <- nrow(ensemblGenes)

Sys.time()
foreach(i=1:numGenes)%dopar%runCisEQTL(i=i,ensemblGenes=ensemblGenes,
                                       exprLocation="data/expr/normalisedCounts/genic/geneExons/byGene_snps1Mb/",
                                       snpLocation="/home/guelfi/eQTL/snps/byGene/",
                                       outputFolder="data/results/genic/geneExons",
                                       genotypeFile="/home/guelfi/plinkOutput/eigenvec",
                                       fullResults="data/results/genic/geneExons/fullResults/")

Sys.time()
stopCluster(cl)
rm(cl)



