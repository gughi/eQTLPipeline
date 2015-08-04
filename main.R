## Author: Manuel Sebastian Guelfi
## email: m.guelfi@ucl.ac.uk
## This is the main script for the eQTL analysis 

## Combined raw counts from genic and intergenic to generate the PEER axes Steps 1 to 
  



## 1. load the data 
  ## 1.1 load the genic data
    exprGenic <- read.csv("data/expr/rawCounts/genic/exprExonIntr.csv", row.names=1)
  ## 1.2 load the intergenic* data
  ## Intergenic data was counts were generated using DerFinder
    load("data/expr/rawCounts/intergenic/allChromosomes.rda")
  ## load data has 2 different objects: annotation,coverage
  ## load the sample info
    load("data/general/sampleInfo.rda")
    
## 2. Filtering
  ## 2.1 Filtering for the genic data
    ## genes that have less than 20 reads in less than 80% of the samples are filtered
    exprGenic <- t(exprGenic)
    ## convert the genes that have NAs
    exprGenic[is.na(exprGenic)]=0
    ## remove genes that not expressed in any gene
    exprGenic <- exprGenic[rowSums(exprGenic>0)>0,]
    ## exprGenic <- exprGenic[rowSums(exprGenic>=20)>(ncol(exprGenic)-((ncol(exprGenic)*20)/100)),]
  ## 2.2 Filtering intergenic data
    ## create matrix we need with the 4 firsts columns that indicate the position of the region
    ## chromosome - start - end - width
    
    exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
    colnames(exprIntergenic) <- c("chr","start","end","width",colnames(coverage))
    rm(annotation,coverage)    
    
    ## select regions with length => 100bp
    exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
    # We don't do filtering since it doesn't make any since for DERFINDER; 
    # regions are detected based on the expression
    

    # load the GC content genic
    GCcontent <- read.delim("data/general/GCcontentExoInt",row.names=1)
    colnames(GCcontent) <- "GCcontent"
    GCcontent$GCcontent <- (GCcontent$GCcontent/100)
    GCcontentTmp <- GCcontent[as.character(rownames(exprGenic)),]
    GCcontentTmp <- as.data.frame(GCcontentTmp)
    rownames(GCcontentTmp) <- rownames(exprGenic)
    colnames(GCcontentTmp) <- "GCcontent"
    rm(GCcontent)
    
    # load the GC content intergenic that has been calculated with BED tools
    GCcontent <- read.delim(pipe("cut -f4,6 data/general/intergenicGCcontent"))
    colnames(GCcontent) <- c("ID","GCcontent")
    rownames(GCcontent) <- GCcontent$ID
    GCcontent$ID <- NULL
    head(GCcontent)
    
    
    ## load the library size
    librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
    
    ## load the gene length for genic counts
    geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
    length <- (geneLength$Gene.End..bp.-geneLength$Gene.Start..bp.)
    geneLength <- cbind(geneLength,length)
    rm(length)
    length <- geneLength[as.character(rownames(GCcontentTmp)),4]    
    
    GCcontentTmp <- rbind(GCcontentTmp,GCcontent)    
    rm(GCcontent)
    
    
    
