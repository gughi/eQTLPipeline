## Author: Manuel Sebastian Guelfi
## email: m.guelfi@ucl.ac.uk
## This is the main script for the eQTL analysis 


## Combined raw counts from genic and intergenic to generate the PEER axes Steps 1 to   

## prerequisites: load the generic function created for the eQTL pipeline
sys.source("genFun.R",
           attach(NULL, name="myenv"))

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
    # creation of identifiers for the regions
    IDs <- paste0("DER",c(1:nrow(exprIntergenic)))
    rownames(exprIntergenic) <- IDs
    rm(IDs)
    
    
    # load the GC content genic and gene length

    geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
    geneLength <- geneLength[as.character(rownames(exprGenic)),c(3,1:2)]
    GCcontent <- GCcalculation(geneLength,genRef="/home/seb/reference/genome37.72.fa")
    length <- (geneLength$Gene.End..bp.-geneLength$Gene.Start..bp.)
    names(length) <- rownames(geneLength)
    GCcontentGenic <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
    colnames(GCcontentGenic)[2] <- "length"
    rm(length,GCcontent)    

    # load the GC content intergenic and region length
    GCcontent <- GCcalculation(exprIntergenic[,1:4],"/home/seb/reference/genome37.72.fa")
    GCcontent <- as.data.frame(cbind(GCcontent[as.character(rownames(exprIntergenic)),],
                      exprIntergenic$width))    
    rownames(GCcontent) <- rownames(exprIntergenic)
    colnames(GCcontent) <- c("GCcontent","length")
    
    ## combined the GCcontent for both quantifications

    GCcontent <- rbind(GCcontentGenic,GCcontent)
    rm(GCcontentGenic)

    ## load the library size
    librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
    
    # combine expression
    expr <- rbind(exprGenic[,as.character(sampleInfo$A.CEL_file)],
                  exprIntergenic[,as.character(sampleInfo$A.CEL_file)])

    # apply cqn
    library(cqn)
    library(scales)

    my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent, sizeFactors = librarySize[colnames(expr),1] , verbose = TRUE)
    
  
    par(mfrow=c(1,2))
    cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
    cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
    RPKM.cqn <- my.cqn$y + my.cqn$offset

    # filter based on the RPKM expression
    # threshold 0.1 that is log2(0.1) because the RPKM are transform in log2
    # we filter genes that have less than 0.1 in 80% of the samples

    
    RPKM.cqn=RPKM.cqn[rowSums(RPKM.cqn>=log2(0.1))>(ncol(RPKM.cqn)-((ncol(RPKM.cqn)*20)/100)),]
    
    # length(grep("ENS",rownames(RPKM.cqn)))
    # 24873 genic regions kept
    # length(grep("DER",rownames(RPKM.cqn)))
    # 28237 intergenic regiond kept

    save(RPKM.cqn,file="data/expr/normalisedCounts/genicIntergenic.rda",compress="bzip2")

        
