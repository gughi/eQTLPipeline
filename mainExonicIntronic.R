## main Exonic+Introns

    setwd("/home/guelfi/eQTLPipeline")
    sink("logExonicIntornic.log")
    nCores <- 15
    cat(paste("Number of cores",nCores,"\n"))
    library(devtools)
    load_all()
    
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
    cat("Processing PUTM \n")
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
    
    cat("load the GC content genic and gene length \n")
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
    cat(paste("number of genes in the analysis for PUTM",length(genesList),"\n"))
    
    cat("Calculating the GC content \n")
    geneLength <- geneLength[as.character(genesList),]
    GCcontent <- GCcalculation(geneLength
                               ,genRef="/home/ukbec/bowtie2Index/genome37.72.fa"
                               ,pathBedtools = "/apps/BEDTools/2.24.0/bin/bedtools")
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
    
    cat("Conditional quantile normalisation \n")
    png(paste0("plots/exonicIntronic/CQNPUTM.jpeg"), type="cairo")
    par(mfrow=c(1,2))
    cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
    cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
    dev.off()
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
    lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
    rownames(lanes) <- gsub("CEL","",rownames(lanes))
    covs <- cbind(covs,librarySize[as.character(rownames(covs))])
    covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
    colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                        "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")
    
    save(RPKM.cqn,PUTM,covs,file="data/expr/normalisedCounts/genic/ExonIntrons/RPKM.cqn.PUTM")
    
    
    rm(list=ls())
    
    load("data/expr/normalisedCounts/genic/ExonIntrons/RPKM.cqn.PUTM")
    
    ##doSwamp(RPKM.cqn,covs)
    
    PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
    PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]
    
    resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                                   "data/expr/normalisedCounts/genic/ExonIntrons/resids.PUTM.rda")
    
    ##doSwamp(resids,covs)
    
    
##############
###  SNIG  ###
##############

    # delete all the objects
    rm(list=ls()) 
    
    nCores <- 15
    load_all()
    
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
    cat("Processing SNIG \n")
    SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),]
    
    # now we select the expression for the SNIG only samples
    expr <- exprGenic[,as.character(SNIG$A.CEL_file)]
    rm(exprGenic)
    
    
    # load the library size
    librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
    librarySize <- librarySize[as.character(SNIG$A.CEL_file),]
    names(librarySize) <- as.character(SNIG$A.CEL_file)
    
    # convert in RPKM
    library(easyRNASeq)
    
    cat("load the GC content genic and gene length \n")
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
    cat(paste("number of genes in the analysis for SNIG",length(genesList),"\n"))
    
    cat("Calculating the GC content \n")
    geneLength <- geneLength[as.character(genesList),]
    GCcontent <- GCcalculation(geneLength
                               ,genRef="/home/ukbec/bowtie2Index/genome37.72.fa"
                               ,pathBedtools = "/apps/BEDTools/2.24.0/bin/bedtools")
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
    
    cat("Conditional quantile normalisation \n")
    png(paste0("plots/exonicIntronic/CQNSNIG.jpeg"), type="cairo")
    par(mfrow=c(1,2))
    cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
    cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
    dev.off()
    RPKM.cqn <- my.cqn$y + my.cqn$offset
    
    SNIG$U.Region_simplified <- NULL
    covs <- SNIG 
    rownames(covs) <- covs$A.CEL_file
    #convert the female and male info in numeric
    covs[covs=="M"]=0
    covs[covs=="F"]=1
    covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
    covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
    covs[,5] <- as.numeric(covs[,5])
    covs[,6] <- as.numeric(covs[,6])
    lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
    rownames(lanes) <- gsub("CEL","",rownames(lanes))
    covs <- cbind(covs,librarySize[as.character(rownames(covs))])
    covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
    colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                        "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")
    
    save(RPKM.cqn,SNIG,covs,file="data/expr/normalisedCounts/genic/ExonIntrons/RPKM.cqn.SNIG")
    
    
    rm(list=ls())
    
    load("data/expr/normalisedCounts/genic/ExonIntrons/RPKM.cqn.SNIG")
    
    ##doSwamp(RPKM.cqn,covs)
    
    PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
    PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]
    
    resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                                   "data/expr/normalisedCounts/genic/ExonIntrons/resids.SNIG.rda")
    
    ##doSwamp(resids,covs)
    
    # delete all the objects
    rm(list=ls())
    
    nCores <- 15
    
    setwd("/home/guelfi/eQTLPipeline")
    cat("starting the split by expression \n")
    writeSH(nameSH="splitByGene.sh",logName="splitByGene",
            cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/parSplitByGeneExonicIntronic.R"),numThreads=(nCores+1))
    
    ### send qsub comand
    system("qsub splitByGene.sh")
    
    ### Run the eQTL analysis
    writeSH(nameSH="runCisEQTL.sh",logName="runCisEQTL",
            cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/parRunCiseQTLExonicIntronic.R"),numThreads=(nCores+1))
    
    ### send qsub comand
    system("qsub runCisEQTL.sh")
    
    writeSH(nameSH="LDsentinalisation.sh",logName="LDsentinalisation",
            cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/sentiExonicIntronic.R"),numThreads=16)
    
    system("qsub LDsentinalisation.sh")
    
    system("perl getAllEQTLsent.pl data/results/genic/ExonIntrons/resMatrixEQTL/sentinalised/PUTM data/results/finaleQTLs/exonicIntronic.PUTM.txt")
    ## Collect for SNIG
    system("perl getAllEQTLsent.pl data/results/genic/ExonIntrons/resMatrixEQTL/sentinalised/SNIG data/results/finaleQTLs/exonicIntronic.SNIG.txt")
    
    # number of eTQL in PUTM at <10% FDR
    system("wc -l data/results/finaleQTLs/exonicIntronic.PUTM.txt")
    # 2290
    
    # number of eTQL in SNIG at <10% FDR
    system("wc -l data/results/finaleQTLs/exonicIntronic.SNIG.txt")
    # 1446

    ## secondary signals SNIG at <1% FDR
    system("cut -f8 -d' ' data/results/finaleQTLs/exonicIntronic.SNIG.txt | sort | grep 2 | wc -l")
    ## 36
    ## tetrary signals SNIG at <1% FDR
    system("cut -f8 -d' ' data/results/finaleQTLs/exonicIntronic.SNIG.txt | sort | grep 3 | wc -l")
    ## 0
    
    ## secondary signals PUTM
    system("cut -f8 -d' ' data/results/finaleQTLs/exonicIntronic.PUTM.txt | sort | grep 2 | wc -l")
    ## 70
    ## tetrary signals PUTM
    system("cut -f8 -d' ' data/results/finaleQTLs/exonicIntronic.PUTM.txt | sort | grep 3 | wc -l")
    ## 2
    
    eQTLPUTM <- read.delim("data/results/finaleQTLs/exonicIntronic.PUTM.txt",sep=" ")
    table(eQTLPUTM$myFDR<0.05)
#     FALSE  TRUE
#     494  1795
    table(eQTLPUTM$myFDR<0.01)
#     FALSE  TRUE
#      1222  1067
    
    
    eQTLSNIG <- read.delim("data/results/finaleQTLs/exonicIntronic.SNIG.txt",sep=" ")
    table(eQTLSNIG$myFDR<0.05)
#     FALSE  TRUE
#     430  1015
    table(eQTLSNIG$myFDR<0.01)
#     FALSE  TRUE
#     886   559
    

    
    
    
    ###########################
    ####    Annotation ########
    ###########################
    
    ### PUTM ###
    
    
    library("biomaRt")
    ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                       dataset="hsapiens_gene_ensembl")
    
    
    geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","gene_biotype","description"),
                       verbose = T,
                       filters="ensembl_gene_id",
                       values=eQTLPUTM$gene, mart=ensembl)
    
    
    
    rownames(geneNames)<- geneNames$ensembl_gene_id
    geneNames$ensembl_gene_id <- NULL
    
    
    eQTLPUTM <- cbind(eQTLPUTM,geneNames[as.character(eQTLPUTM$gene),])
    
    save(eQTLPUTM,file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
    write.csv(eQTLPUTM,file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.csv",row.names=F)    
    
    ### SNIG ###
    
    rm(geneNames)
    
    geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","gene_biotype","description"),
                       verbose = T,
                       filters="ensembl_gene_id",
                       values=eQTLSNIG$gene, mart=ensembl)
    
    
    
    rownames(geneNames)<- geneNames$ensembl_gene_id
    geneNames$ensembl_gene_id <- NULL
    
    
    eQTLSNIG <- cbind(eQTLSNIG,geneNames[as.character(eQTLSNIG$gene),])
    
    save(eQTLSNIG,file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
    write.csv(eQTLSNIG,file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.csv",row.names=F)    
    
    sink()
    
    
    #########################
    ### Annotation by SNP ###
    #########################
    
    rm(list=ls())
    
    load(file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
    eQTLPUTM <-  annSNP(eQTLPUTM)
    head(eQTLPUTM)
    
    save(eQTLPUTM,file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
    write.csv(eQTLPUTM,file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.csv",row.names=F)    
    
    load(file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
    eQTLSNIG <-  annSNP(eQTLSNIG)
    
    save(eQTLSNIG,file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
    write.csv(eQTLSNIG,file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.csv",row.names=F)    
    
    
    ################################################ 
    #### check overlapping genes between tissues ###
    ################################################
    
    
    load(file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
    load(file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
    
    library(gplots)
    ## FDR10%
    venn(list(PUTM=unique(eQTLPUTM$gene),SNIG=unique(eQTLSNIG$gene)))
    
    ## FDR5%    
    eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.05),]
    eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.05),]
    
    venn(list(PUTM=unique(eQTLPUTM$gene),SNIG=unique(eQTLSNIG$gene)))
    
    
    eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.01),]
    eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]
    
    venn(list(PUTM=unique(eQTLPUTM$gene),SNIG=unique(eQTLSNIG$gene)))
    
    
    
    
    