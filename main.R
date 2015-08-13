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
    exprGenic <- exprGenic[rowSums(exprGenic>0)>0,as.character(sampleInfo$A.CEL_file)]
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

    # load the library size
    librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
    librarySize <- librarySize[as.character(sampleInfo$A.CEL_file),]
    names(librarySize) <- as.character(sampleInfo$A.CEL_file)

    # convert in RPKM
    library(easyRNASeq)

    # load the GC content genic and gene length

    geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
    geneLength <- geneLength[as.character(rownames(exprGenic)),c(3,1:2)]
    length <- (geneLength$Gene.End..bp.-geneLength$Gene.Start..bp.)
    names(length) <-  as.character(rownames(exprGenic))
    stopifnot(identical(colnames(exprGenic),names(librarySize)))
    stopifnot(identical(rownames(exprGenic),names(length)))              
    RPKM.std <- RPKM(exprGenic, NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

    ## filtering
    RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
    genesList <- rownames(RPKM.std) 
    rm(RPKM.std,length)

    geneLength <- geneLength[as.character(genesList),]
    GCcontent <- GCcalculation(geneLength,genRef="/home/seb/reference/genome37.72.fa")
    length <- (geneLength$Gene.End..bp.-geneLength$Gene.Start..bp.)
    names(length) <- rownames(geneLength)
    GCcontentGenic <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
    colnames(GCcontentGenic)[2] <- "length"
    exprGenic <- exprGenic[genesList,]
    rm(length,GCcontent,genesList,geneLength)

    
    

    # load the GC content intergenic and region length
        
    length <- exprIntergenic$width    
    names(length) <- rownames(exprIntergenic)
    stopifnot(identical(colnames(as.matrix(exprIntergenic[,as.character(sampleInfo$A.CEL_file)]))
                        ,names(librarySize)))
    stopifnot(identical(rownames(exprIntergenic),names(length)))              
                                  
    #convert in RPKM
    RPKM.std <- RPKM(as.matrix(exprIntergenic[,as.character(sampleInfo$A.CEL_file)])
                     , NULL, 
                     lib.size=librarySize, 
                     feature.size=length)

    RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
    genesList <- rownames(RPKM.std) 
    rm(RPKM.std,length)
    exprIntergenic <- exprIntergenic[as.character(genesList),]
  
    GCcontent <- GCcalculation(exprIntergenic[,1:4],"/home/seb/reference/genome37.72.fa")
    GCcontent <- as.data.frame(cbind(GCcontent[as.character(rownames(exprIntergenic)),],
                      exprIntergenic$width))    
    rownames(GCcontent) <- rownames(exprIntergenic)
    colnames(GCcontent) <- c("GCcontent","length")
    
    ## combined the GCcontent for both quantifications

    GCcontent <- rbind(GCcontentGenic,GCcontent)
    rm(GCcontentGenic)

    # combine expression
    expr <- rbind(exprGenic[,as.character(sampleInfo$A.CEL_file)],
                  exprIntergenic[,as.character(sampleInfo$A.CEL_file)])

    # apply cqn
    library(cqn)
    library(scales)

    stopifnot(identical(colnames(expr),names(librarySize)))
    stopifnot(identical(rownames(expr),rownames(GCcontent)))              

    my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent, sizeFactors = librarySize , verbose = TRUE)
    
  
    par(mfrow=c(1,2))
    cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
    cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
    RPKM.cqn <- my.cqn$y + my.cqn$offset

    # filter based on the RPKM expression
    # threshold 0.1 that is log2(0.1) because the RPKM are transform in log2
    # we filter genes that have less than 0.1 in 80% of the samples


    # length(grep("ENS",rownames(RPKM.cqn)))
    # 25985 genic regions kept
    # length(grep("DER",rownames(RPKM.cqn)))
    # 42577 intergenic regiond kept
    
    # save results
    save(RPKM.cqn,file="data/expr/normalisedCounts/genicIntergenic.rda",compress="bzip2")

    
    # PCA PC1 vs PC2 
    PCAres<- prcomp(t(RPKM.cqn))
    par(mfrow=c(1,1))
    PUTM <- sampleInfo[sampleInfo[, 6] == "PUTM",]
    SNIG <- sampleInfo[sampleInfo[, 6] == "SNIG",]
    
    plot(PCAres, main="PCA axis simple quantification (PUTM + SNIG)")
    plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue",xlab="PC1",ylab="PC2" )
    points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
    points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
    legend("bottomright", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    
    
    rm(PUTM,SNIG,PCAres)
        
    # heatmap to check correlation between PCs adn known factor (e.g. sex, age, etc...)
    
    # This is the code used to produced the covs
     covs <- sampleInfo
     rownames(covs) <- covs$A.CEL_file
     #convert the female and male info in numeric
     covs[covs=="M"]=0
     covs[covs=="F"]=1
     covs[covs=="PUTM"]=1
     covs[covs=="SNIG"]=2
     covs <- as.data.frame(apply(covs[,c(2:6,8:10)], 2, as.factor))
     covs[,c(1:5,8)] <- as.data.frame(apply(covs[,c(1:5,8)], 2, as.numeric))
     covs[,6] <- as.numeric(covs[,6])
     covs[,7] <- as.numeric(covs[,7])
     lanes <- read.csv("/home/seb/expressionData/QCmetrics.csv",row.names=8)
     rownames(lanes) <- gsub("CEL","",rownames(lanes))
     covs <- cbind(covs,librarySize[as.character(rownames(covs))])
     covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
     colnames(covs) <- c("Age","PMI","RIN","Gender","Region","CODE","OVation_Batch",
                        "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")
  

    doSwamp(RPKM.cqn=RPKM.cqn,covs)


    save(RPKM.cqn,covs,file="data/general/RPKMCQNcovs.rda")    

    
    # PEER look at the script for the optimisation    
    
    # How I select the best PEER I basically calculate the rsquared for all the PEER test,
    # from each test calculate the maximum and then I select the PEER test that has minimum rsquared beetween all the tests  
    load("testPEER/correlation.rda")
    print(paste("Test with minimun correlation",match(min(apply((corPEER^2),2,max)),apply((corPEER^2),2,max))))
    
    load("data/general/RPKMCQNcovs.rda")
    PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
    correPlot(PEERRNDPEER18,covs[rownames(PEERRNDPEER18),],"Correlation 'best' PEER and known factors")


    # delete all the objects
    rm(list=ls())
    ## Now I separated the analysis for each quantification
    
    ## Now we correct for PEER using simple quantification Exons+Introns
    exprGenic <- read.csv("data/expr/rawCounts/genic/exprExonIntr.csv", row.names=1)
    # load the sample info to get the IDs for each tissue
    load("data/general/sampleInfo.rda")
  
    # transpose the expression matrix this to have the format of: In the rows the observations(genes) and in the columns the samples
    exprGenic <- t(exprGenic)
    ## convert the genes that have NAs
    exprGenic[is.na(exprGenic)]=0
    ## remove genes that not expressed in any gene
    exprGenic <- exprGenic[rowSums(exprGenic>0)>0,]

    PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
    
    # now we select the expression for the PUTM only samples
    expr <- exprGenic[,as.character(PUTM$A.CEL_file)]
    rm(exprGenic)

    
    # load the library size
    librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
    librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
    names(librarySize) <- as.character(PUTM$A.CEL_file)
    
    # convert in RPKM
    library(easyRNASeq)
    
    # load the GC content genic and gene length
    
    geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
    geneLength <- geneLength[as.character(rownames(expr)),c(3,1:2)]
    length <- (geneLength$Gene.End..bp.-geneLength$Gene.Start..bp.)
    names(length) <-  as.character(rownames(expr))
    stopifnot(identical(colnames(expr),names(librarySize)))
    stopifnot(identical(rownames(expr),names(length)))              
    
    RPKM.std <- RPKM(expr, NULL, 
                     lib.size=librarySize, 
                     feature.size=length)
    
    ## filtering
    RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
    genesList <- rownames(RPKM.std) 
    rm(RPKM.std,length)


    # now we calculate the GC content 
    geneLength <- geneLength[as.character(genesList),]
    GCcontent <- GCcalculation(geneLength
                               ,genRef="/home/seb/reference/genome37.72.fa")
    length <- (geneLength$Gene.End..bp.-geneLength$Gene.Start..bp.)
    names(length) <- rownames(geneLength)
    GCcontent <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
    colnames(GCcontent)[2] <- "length"
    # we update gene expression with the filtered genes
    expr <- expr[genesList,]
    rm(length,geneLength,genesList)

    ## load the library size
    
    library(cqn)
    library(scales)
    
    ## CQN
    stopifnot(identical(colnames(expr),names(librarySize)))
    stopifnot(identical(rownames(expr),rownames(GCcontent)))              

    my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent,sizeFactors=librarySize, verbose = TRUE)

    par(mfrow=c(1,2))
    cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
    cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
    RPKM.cqn <- my.cqn$y + my.cqn$offset
    
    PUTM$U.Region_simplified <- NULL
    covs <- PUTM 
    rownames(covs) <- covs$A.CEL_file
    #convert the female and male info in numeric
    covs[covs=="M"]=0
    covs[covs=="F"]=1
    covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
    covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
    covs[,5] <- as.numeric(covs[,5])
    covs[,6] <- as.numeric(covs[,6])
    lanes <- read.csv("/home/seb/expressionData/QCmetrics.csv",row.names=8)
    rownames(lanes) <- gsub("CEL","",rownames(lanes))
    covs <- cbind(covs,librarySize[as.character(rownames(covs))])
    covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
    colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                    "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")

    save(RPKM.cqn,PUTM,covs,file="data/expr/normalisedCounts/genic/Exon+Introns/RPKM.cqn.PUTM")
    
    ### Residual correction ###

    rm(list=ls())

    load("data/expr/normalisedCounts/genic/Exon+Introns/RPKM.cqn.PUTM")
  
    doSwamp(RPKM.cqn,covs)

    PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
    PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

    resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                         "/home/seb/projectsR/eQTLPipeline/data/expr/normalisedCounts/genic/Exon+Introns/resids.PUTM.rda")
    
    doSwamp(resids,covs)

    



        