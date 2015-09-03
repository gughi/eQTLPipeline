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
clusterExport(cl,c("LDsentinalisation","read.table.rows"))
registerDoParallel(cl)
getDoParWorkers()
## path we put the results
pathFinalSentinalised <-"data/results/genic/geneExons/resMatrixEQTL/sentinalised/FDR10/"
## path where get the unsentinalised eQTLs
pathUnsentinalised <- "data/results/genic/geneExons/resMatrixEQTL/"



ensemblGenes <- read.delim(pipe(paste0("ls ",pathUnsentinalised,
                                       "PUTM/ > ",pathUnsentinalised,"genes.tmp.txt; ls ",pathUnsentinalised,
                                       "SNIG/ >> ",pathUnsentinalised,"genes.tmp.txt; sort ",pathUnsentinalised,
                                       "genes.tmp.txt | uniq; rm ",pathUnsentinalised,
                                       "genes.tmp.txt")),header=F)


Sys.time()
foreach(i=1:nrow(ensemblGenes))%dopar%LDsentinalisation(i=i,
                                                        ensemblGenes=ensemblGenes,
                                                        pathFinalSentinalised=pathFinalSentinalised,
                                                        pathUnsentinalised=pathUnsentinalised,
                                                        FDRthr=0.10,
                                                        exprLocation="data/expr/normalisedCounts/genic/geneExons/byGene_snps1Mb/",
                                                        snpLocation="/home/guelfi/eQTL/snps/byGene/",
                                                        genotypeFile="/home/guelfi/plinkOutput/eigenvec",
                                                        tmpFolder="tmp/")
##foreach(i=1:20)%dopar%splitExprByGene(i,ensemblRef,ensemblGenes,PUTM,exprPUTM,SNIG,exprSNIG)
Sys.time()
stopCluster(cl)

