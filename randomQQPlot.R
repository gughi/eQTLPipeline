### This script select 100 random pvalues for each gene for the eQTL results and then plot the QQplot 

setwd("/home/guelfi/eQTLPipeline/")

getRandomPval <- function(i,path,genes,numpvals=100){
  ## number of genes
  eQTLs <- read.delim(paste0(path,genes[i,1]))
  rand <- 0
  if(nrow(eQTLs) > numpvals)
  {
    rand <- sample(2:nrow(eQTLs), numpvals, replace=FALSE)
  }
  else 
  {
    rand <- 1:nrow(eQTLs)
  }
  
  return(eQTLs[rand,5])
} 

library(devtools)
library(doParallel)
library(foreach)
load_all()

cl <- makeCluster(7)
registerDoParallel(cl)
getDoParWorkers()

source(file="/gpfs/Caprica/home/seb/Rscripts/mikeFunctions.R")

path <- "data/results/genic/ExonIntrons/fullResults/PUTM/"
genes <- read.delim(pipe(paste0("ls ",path)),header=F)

Sys.time()
pvalues <-foreach(i=1:nrow(genes),.combine=c)%dopar%getRandomPval(i,path,genes)
Sys.time()
stopCluster(cl)
rm(cl)
save(pvalues,file="tmp/pValGeneExoIntron.PUTM.rda")
##load(file="tmp/pValGeneExo.PUTM.rda")

png(paste0("plots/qqplotExonicIntronicPUTM.jpeg"), type="cairo")
qq.plot(pvalues)
dev.off()

## intergenic

library(devtools)
library(doParallel)
library(foreach)
load_all()

cl <- makeCluster(4)
registerDoParallel(cl)
getDoParWorkers()

source(file="/gpfs/Caprica/home/seb/Rscripts/mikeFunctions.R")

path <- "data/results/intergenic/fullResults/SNIG/"
genes <- read.delim(pipe(paste0("ls ",path)),header=F)

Sys.time()
pvalues <-foreach(i=1:nrow(genes),.combine=c)%dopar%getRandomPval(i,path,genes,10)
Sys.time()
stopCluster(cl)
rm(cl)
save(pvalues,file="tmp/pValIntergenic.SNIG.rda")

png(paste0("plots/qqplotIntergenicSNIG.jpeg"), type="cairo")
qq.plot(pvalues)
dev.off()


