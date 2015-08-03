## Author: Manuel Sebastian Guelfi
## email: m.guelfi@ucl.ac.uk
### This is the main script for the eQTL analysis 

## Combined raw counts from genic and intergenic to generate the PEER axes Steps 1 to 
  



## 1. load the data 
  ## 1.1 load the genic data
    exprGenic <- read.csv("expr/rawCounts/genic/exprExonIntr.csv", row.names=1)
  ## 1.2 load the intergenic* data
  ## Intergenic data was counts were generated using DerFinder
    load("expr/rawCounts/intergenic/allChromosomes.rda")
  ## load data has 2 different objects: annotation,coverage

## 2. Filtering
  ## 2.1 Filtering for the genic data
    ## genes that have less than 20 reads in less than 80% of the samples are filtered
    exprGenic <- t(exprGenic)
    ## convert the genes that have NAs
    exprGenic[is.na(exprGenic)]=0
    ## remove genes
    exprGenic <- exprGenic[rowSums(exprGenic>=20)>(ncol(exprGenic)-((ncol(exprGenic)*20)/100)),]
  ## 2.2 Filtering intergenic data
    ## create matrix we need with the 4 firsts columns that indicate the position of the region
    ## chromosome - start - end - width
    
    exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
    colnames(expr) <- c("chr","start","end","width",colnames(coverage))
    rm(annotation,coverage)    
    
    ## select regions with length => 100bp
    exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
    ## genes that have less than 20 reads in less than 80% of the samples are filtered  
    exprIntergenic <- exprIntergenic[rowSums(exprIntergenic>=20)>(ncol(exprIntergenic)-((ncol(exprIntergenic)*20)/100)),]
    
    





