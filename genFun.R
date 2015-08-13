## This script contains the general gunction needed for the eQTL analysis

#' Calculate the GC content
#' 
#' @param data frame with the regions with the following fixed structure: first column is the chr, second column is the start position, third column is the stop position, forth column is the Identifier
#' @param the location of the .fa file of the genome reference to calculate the GC content
#' @note this function need bedtools installed in your system
#' @return data frame with one column as GC content, the rownames are the identifiers
GCcalculation <- function (region,genRef) {
  ##We now define the BED file to then use it to calculate GC content
  print("creting the BED file")
  tmpBED <- tempfile("GCcont", fileext = ".BED")
  BED <- paste(gsub("chr","",region[,1]),region[,2],region[,3],rownames(region), sep="\t")
  write.table(data.frame(BED), file = tmpBED, row.names = F, 
              col.names = F, quote = F)
  rm(BED)
  
  tmpGCcon <- tempfile("GCcont")  
  cmd <- paste0("bedtools nuc -fi ",genRef," -bed ",tmpBED," > ",tmpGCcon)
  
  print("executing the bedtools")
  system(cmd)
  rm(cmd)
  
  print("collecting the GC content")
  GCcontent <- read.delim(pipe(paste("cut -f4,6", tmpGCcon)))
  colnames(GCcontent) <- c("ID","GCcontent")
  rownames(GCcontent) <- GCcontent$ID
  GCcontent$ID <- NULL 
  GCcontent
}


doSwamp <- function(RPKM.cqn,covs)
{
  library(swamp)
  RPKM.cqn.tmp <- t(RPKM.cqn)
  expr.data.o <- RPKM.cqn.tmp[order(rownames(RPKM.cqn.tmp)),]
  
  ## Order your data in the same way
  traits.o <- covs[as.character(rownames(expr.data.o)),]
  rownames(traits.o) <- rownames(expr.data.o)
  
  
  print(head(covs))
  stopifnot(identical(rownames(traits.o),rownames(expr.data.o)))
  
  
  #Do PCA analysis with up to 10 axis
  res1 <- prince(t(expr.data.o),traits.o,top=15,permute=TRUE)
  #PLOT WHAT YOU NEED
  prince.plot(prince=res1)
  
  #Plot the variation included by each PCA axis
  res2 <- prince.var.plot(t(expr.data.o),show.top=50,npermute=10)
  
  hca.plot(as.matrix(t(expr.data.o)),as.data.frame(traits.o))

}

#' Function that calculates hidden factors in the experiment usin PEER 
#' @param RPKM.cqn = The expression normalised
#' @param nFactors the number of hidden factor we want to discover
#' @param covs the covariates to exclude from PEER analysis
#' @param outputFile where the results are going to be saved
doPEER <- function(RPKM.cqn,nFactors,covs,outputFile)
{
  
  library(peer,lib="/home/seb/Rlibrary/")
  
  
  ## WE DON'T INCLUDE THE GENETIC COVARIATS
  ## geneticPCA <- read.table.rows("/home/guelfi/plinkOutput/eigenvec", keepRows=covs$U.SD_No, sep=" ",header=F)
  ## geneticPCA <- geneticPCA[covs$U.SD_No,]
  ## rownames(geneticPCA) <- rownames(covs)
  
  ## WE INCLUDE AGE,GENDER AND TISSUE
  covs <- covs[as.vector(rownames(RPKM.cqn)),]
  stopifnot(identical(rownames(covs),rownames(RPKM.cqn)))
  
  ## create the PEER model
  model=PEER()
  ## set to obtain only the first 15 hidden factors
  PEER_setNk(model, nFactors)
  ## set the covariants, step 4
  PEER_setCovariates(model, data.matrix(covs))
  PEER_setPhenoMean(model,as.matrix(RPKM.cqn))
  
  
  ## check if the covariants have been updatef properly
  print(paste("dimension of the expression matrix",print(dim(PEER_getPhenoMean(model)))))
  print(paste("dimension of the covariate matrix",print(dim(PEER_getCovariates(model)))))
  ##head(PEER_getPhenoMean(model))
  
  ## run PEER, step 5
  PEER_update(model)
  
  ## get the hidden factors
  factor <- PEER_getX(model)
  head(factor)
  rownames(factor) <- rownames(RPKM.cqn)
  
  write.csv(factor,file=outputFile)
}


correPlot <- function (mat1,mat2,xlab)
{
  library(gplots)
  linp<-matrix(ncol=ncol(mat1),nrow=ncol(mat2))  
  rownames(linp)<-colnames(mat2)
  colnames(linp)<-colnames(mat1)
  rsquared<-matrix(ncol=ncol(mat1),nrow=ncol(mat2))
  rownames(rsquared)<-colnames(mat2)
  colnames(rsquared)<-colnames(mat1)
  for (i in 1:ncol(mat2)){
    for (j in 1:ncol(mat1)){
      fit<-lm(mat1[,j]~mat2[,i])
      s<-summary(fit)
      linp[i,j]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
      rsquared[i,j]<-s$r.squared[1]
    }}
  
  
  smallest=-20
  linp10<-log10(linp)
  linp10<-replace(linp10,linp10<=smallest,smallest)
  print(linp10)
  
  rsquaredTMP <- rsquared
  rsquaredTMP[which(linp10>-5)] = NA
  rsquaredTMP <-round(rsquaredTMP,digits=2)
  heatmap.2(linp10,Colv=F,Rowv=F,dendrogram="none",trace="none",symbreaks=F,symkey=F,breaks=seq(-20,0,length.out=100),key=T,col=heat.colors(99),
            cexRow=1,cexCol=1,colsep=NULL,rowsep=NULL,sepcolor=sepcolor,sepwidth=sepwidth,
            main="",labCol=paste(1:ncol(linp10),sep=""),margins=c(5,7),labRow=,xlab=xlab,
            cellnote=rsquaredTMP,notecol="black",notecex=1) 
#   heatmap.2(linp10,Colv=F,Rowv=F,dendrogram="none",trace="none",symbreaks=F,symkey=F,breaks=seq(-20,0,length.out=100),key=T,col=heat.colors(99),
#             cexRow=1,cexCol=1,colsep=NULL,rowsep=NULL,sepcolor=sepcolor,sepwidth=sepwidth,
#             main="",labCol=paste(1:ncol(linp10),sep=""),margins=c(5,7),labRow=,xlab=xlab,
#             cellnote=matrix(ncol=ncol(linp),nrow=nrow(linp)),notecol="black",notecex=1) 
#   
  
  
}


#' Function that does that returns the residuals 
#' @param exprThe expression matrix nxm where n are the samples and m the genes
#' @param covariates to use for the correction
#' @param covs the covariates to exclude from PEER analysis
#' @param outputFile where the results are going to be saved 
doResidualCorrection <- function(expr,covs,outputFile)
{
  
  ## load the knowing factors (age and sex ) and unknown factors (PEER axes)
  ## Residual correction with linear model, step 3
  print("Performing the residual correction")
  stopifnot(identical(rownames(expr),rownames(covs)) )
  resids <- apply(expr, 2, function(y){
    lm( y ~ . , data=covs)$residuals 
  })
  
  print("Saving the residuals")
  save(resids,file=paste0(outputFile))
  resids
}




splitExprByGene <- function(i,ensemblRef,ensemblGenes,PUTM,exprPUTM,SNIG,exprSNIG,pathResExpr)
{
  ## general functions ##
  sys.source("/home/adai/scripts/common_functions.R",
             attach(NULL, name="myenv"))
  
  
  ## load the ensembl file with ENSG, ENST, gene start and stop and chr
  ## ensemblRef <- read.delim(file="/home/seb/eQTL/ensemblRef.txt", as.is=T,header=T)
  
  ##load the entire list of genes
  ## ensemblGenes <- read.delim(file="/home/seb/eQTL/ensemblGenes.txt", as.is=T,header=T)
  ensemblRef.gene  <- ensemblGenes[i, ]; rm(i)
  geneID <- ensemblRef.gene$Ensembl.Gene.ID
  ## get the ensemble geneID
  gene.info <- ensemblRef[ which(ensemblRef$Ensembl.Gene.ID == geneID), ]
  
  print(geneID)
  ## "ENSG00000266195"
  
  ## filenames and folders ##
  ##pathResExpr <- "/home/seb/eQTL/expr/simpleQuantificationExonIntrons/CQNINT/13PEER/"
  dir.create(pathResExpr, showWarnings=FALSE )
  ##dir.create("/home/seb/eQTL/expr/simpleQuantificationExonIntrons/CQN/9PEER/chunks_log/", showWarnings=FALSE )
  
  fn.rda <- paste0(pathResExpr,"/byGene_snps1Mb/", geneID, ".rda")
  dir.create( dirname(fn.rda), showWarnings=FALSE )
  
  ## Read in tissues
  tissues <- c("PUTM","SNIG")
  expr <- vector(mode="list", length(tissues))
  names(expr) <- tissues
  
  ## change the sample ID to the individual ID to match the genotyped data
  #samples <- rownames(exprPUTM)
  # indID <- read.table.rows("/home/seb/phenotype/PUTMinfo.csv", keepRows=samples, sep=",")
  # IDs <- indID[samples,6]
  IDs <- PUTM[which(PUTM$A.CEL_file %in% as.character(rownames(exprPUTM))),1]
  IDs <- gsub("/","_",IDs)
  rownames(exprPUTM) <- IDs
  rm(indID,IDs)
#   samples <- rownames(exprSNIG)
#   indID <- read.table.rows("/home/seb/phenotype/SNIGinfo.csv", keepRows=samples, sep=",")
#   IDs <- indID[samples,6]
#   rownames(exprSNIG) <- IDs
  IDs <- SNIG[which(SNIG$A.CEL_file %in% as.character(rownames(exprSNIG))),1]
  IDs <- gsub("/","_",IDs)
  rownames(exprSNIG) <- IDs
  rm(IDs)
  
  
  expr[["SNIG"]] <- as.matrix(exprSNIG[,is.element(colnames(exprSNIG),ensemblRef.gene$Ensembl.Gene.ID)])
  expr[["PUTM"]] <- as.matrix(exprPUTM[,is.element(colnames(exprPUTM),ensemblRef.gene$Ensembl.Gene.ID)])
  
  ## check if the gene is expressed in one tissu or 
  if (length(expr[["SNIG"]])==0 && length(expr[["PUTM"]])==0)
  {
    cat(geneID,"\n",file=paste0(pathResExpr,"/genesNoExpressed"), append=T)
    stop("Gene not expressed")
  }
  for(tissue in tissues ){
    if (length(expr[[tissue]])==0)
    {
      cat(geneID,"\n",file=paste0(pathResExpr,"/genesNoExpressed",tissue), append=T)
      tissues <- tissues[-which(tissues==tissue)]
    }
    ## this case assign the column name to maintain the structure of the matrix, specially to keep the name of the transcript
    
    if (ncol(expr[["SNIG"]])==1)
    {
      colnames(expr[["SNIG"]]) <- intersect(colnames(exprSNIG),ensemblRef.gene$Ensembl.Gene.ID)
    }
    if (ncol(expr[["PUTM"]])==1)
    {
      colnames(expr[["PUTM"]]) <- intersect(colnames(exprPUTM),ensemblRef.gene$Ensembl.Gene.ID)
    }
  }
  
  save(expr, gene.info,
       file=fn.rda, compress="bzip2", ascii=T)
  rm(fn.rda)
  
}


runCisEQTL <- function(i,ensemblGenes,exprLocation,snpLocation,outputFolder,genotypeFile)
{
  
  ## general functions ##
  sys.source("/home/adai/scripts/common_functions.R",
             attach(NULL, name="myenv"))
  
  
  ##load the entire list of genes
  ##ensemblGenes <- read.delim(pipe("ls /home/seb/eQTL/expr/simpleQuantificationExonIntrons/CQN/11PEER/byGene_snps1Mb/"),header=F)
  ##ensemblGenes <- read.delim("/home/guelfi/expr180/simpleQuantificationExonIntrons/CQN/11PEER/geneList",header=F)
  geneID  <- ensemblGenes[i,1 ]; rm(i)
  geneID  <- sub(".rda","",geneID)
  print(geneID)
  
  ## load snps and expression information
  # fn.rda <- paste0("/home/seb/eQTL/expr/simpleQuantificationExonIntrons/CQNINT/13PEER/byGene_snps1Mb/", geneID, ".rda")
  fn.rda <- paste0(exprLocation, geneID, ".rda")
  load(fn.rda)
  #fn.rda <- paste0("/home/seb/eQTL/snps/byGene/", geneID, ".rda")
  fn.rda <- paste0(snpLocation, geneID, ".rda")
  load(fn.rda)
  rm(fn.rda, gene.info, t.map)
  
  ##dirRes <- "/home/seb/eQTL/resEQTLs/simpleQuantificationExonIntrons/CQN/13PEER/FDR10/"
  dir.create(file.path(outputFolder),showWarnings=FALSE)
  dir.create(file.path(outputFolder,"resMatrixEQTL/"),showWarnings=FALSE)
  dir.create(file.path(paste0(outputFolder,"resMatrixEQTL/"), "SNIG"),showWarnings=FALSE)
  dir.create(file.path(paste0(outputFolder,"resMatrixEQTL/"), "PUTM"),showWarnings=FALSE)
  
  
  tissues <- c("PUTM","SNIG")
  
  if(markers=="No polymorphisms found within +/- 1Mb of TSS or TES")
  {
    stop("No polymorphisms found within +/- 1Mb of TSS or TES")
  }
  
  library(MatrixEQTL)
  
  
  for( tissue in tissues ){
    
    ## check if expr and markers have the same number of individuals, basically if all the individuals 
    ## have expression data for that gene but also markers
    print(tissue)
    my.expr0 <- as.matrix( expr[[tissue]] )
    
    my.expr0 <- as.matrix( my.expr0[ , which(colMeans( is.na(my.expr0) ) != 1) ])
    
    ## check if there is no expression for the tissue
    if(ncol(my.expr0)<1){next}
    ### ATTENTION the following command is used because some genotyped data are missing
    my.markers0 <- as.matrix(markers[ , is.element(colnames(markers),rownames(my.expr0))])
    
    
    ### add the genetic covariants
    # my.covTMP <- read.table.rows(paste0("/home/seb/plinkOutput/eigenvec"), keepRows=rownames(my.expr0), sep=" ",header=F)
    my.covTMP <- read.table.rows(paste0(genotypeFile), keepRows=rownames(my.expr0), sep=" ",header=F)
    my.cov0 <- as.matrix(my.covTMP[colnames(my.markers0),2:4])
    my.cov0 <- t(my.cov0)
    ##print(my.cov0)
    
    
    ## transpose needed to have the individuals as columns
    my.expr0 <- as.matrix(my.expr0[colnames(my.markers0),])
    ## in case number of columns is equal to 1 we assign again the name of the column again
    if (ncol(my.expr0)==1)
    {
      colnames(my.expr0) <- colnames(expr[[tissue]])
    }
    
    my.expr0 <- t(my.expr0)
    
    ##print(head(my.expr0))
    
    stopifnot( identical( colnames(my.expr0), colnames(my.markers0) ) )
    stopifnot( identical( colnames(my.cov0), colnames(my.markers0) ) )
    ## create the object SlicedData
    my.expr    <- SlicedData$new()
    my.expr$CreateFromMatrix( my.expr0 )
    my.markers <- SlicedData$new()
    my.markers$CreateFromMatrix(my.markers0)
    my.cov <- SlicedData$new()
    my.cov$CreateFromMatrix(my.cov0)
    rm(my.expr0, my.markers0,my.cov0)
    
    ## outputfile
    outputFile=paste0(outputFolder,"resMatrixEQTL/",tissue,"/",geneID)
    
    
    store <- Matrix_eQTL_main( my.markers, my.expr, my.cov,output_file_name = NULL,pvOutputThreshold=1e-1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
    
    
    ## we calculate the FDR without taking the min as the method used by Shabalin
    pval <- store$all$eqtls$pvalue
    myFDR <- sapply(pval, function(pval) p.adjust(pval, method="fdr", n =  nrow(markers.info))) 
    rm(pval)
    my.eQTLstmp <- cbind(store$all$eqtls, myFDR)
    ##print(head(my.eQTLstmp))    
    ## Using a FDR threshold of 5%
    my.eQTLs <- my.eQTLstmp[which(my.eQTLstmp$myFDR <= 0.10),]
    rm(my.eQTLstmp) 
    ##print(head(my.eQTLs))
    
    print(nrow(my.eQTLs))
    if (nrow(my.eQTLs)>0)
    {
      write.table(my.eQTLs,outputFile,row.names=F)
    }
    
    
    
    
  }
  
  rm(my.expr,my.markers,my.cov, markers.info)
}




LDsentinalisation <- function(i,ensemblGenes,pathFinalSentinalised,pathUnsentinalised,FDRthr,exprLocation,snpLocation,genotypeFile,tmpFolder)
{
  ## general functions ##
  sys.source("/home/adai/scripts/common_functions.R",
             attach(NULL, name="myenv"))
  
  
  geneID  <- ensemblGenes[i,1 ]; rm(i)
  print(geneID)
  
  dir.create(file.path(pathFinalSentinalised),showWarnings=FALSE)
  dir.create(file.path(pathFinalSentinalised, "SNIG"),showWarnings=FALSE)
  dir.create(file.path(pathFinalSentinalised, "PUTM"),showWarnings=FALSE)
  
  ## initialise my.cov0
  
  
  
  fn.rda <- paste0(exprLocation, geneID, ".rda")
  load(fn.rda)
  fn.rda <- paste0(snpLocation, geneID, ".rda")
  load(fn.rda)
  
  rm(fn.rda)
  tissues <- c("SNIG","PUTM")
  
  library(MatrixEQTL)
  for( tissue in tissues ){
    print(tissue)
    
    if(!file.exists(paste0(pathUnsentinalised,tissue,"/",geneID)))
    {
      next
    }
    eQTLs <- read.delim(paste0(pathUnsentinalised,tissue,"/",geneID),sep=" ")
    ## j is the degree subsignals
    j=1
    my.cov0 <- NA
    while(TRUE)
    {
      
      ## check if expr and markers have the same number of individuals, basically if all the individuals have expression data for that gene but also markers
      ## check whether there is any other eQTL
      if (j>1)
      { 
        if(!file.exists(paste0(tmpFolder,tissue,"/",geneID)))
        {
          break
        }
        ##system(paste0("rm /home/guelfi/expr180/tmp/",tissue,"/",geneID))
        eQTLs <- read.delim(paste0(tmpFolder,tissue,"/",geneID),sep=" ")
        
      }
      else
      {
        system(paste("echo","snps","gene","t-stat","pvalue","FDR","beta","myFDR","degree", paste0(">> ",pathFinalSentinalised,tissue,"/",geneID)))
      }
      
      if (min(eQTLs$myFDR)>FDRthr)
      {
        print(min(eQTLs$myFDR))
        break
      }
      
      ##if(eQTLs$snps=="No significant associations were found." | nrow(eQTLs)<1)
      ##{
      if(file.exists(paste0(tmpFolder,tissue,"/",geneID)))
      {
        system(paste0("rm ",tmpFolder,tissue,"/",geneID))
      }
      ##break
      ##}
      
      SNP <- eQTLs$snps[which(eQTLs$myFDR==min(eQTLs$myFDR))[1]]        
      leadEQTL <- eQTLs[which(eQTLs$myFDR==min(eQTLs$myFDR))[1],]
      
      print(eQTLs$snps[which(eQTLs$myFDR==min(eQTLs$myFDR))])
      ##system(paste("echo",leadEQTL$SNP,leadEQTL$gene,leadEQTL$beta,leadEQTL$t.stat,leadEQTL$p.value,leadEQTL$FDR, j, paste0(">> /home/guelfi/expr180/tmp/sentinalised",tissue)))
      tmpRes <- data.frame(cbind(eQTLs[which(eQTLs$myFDR==min(eQTLs$myFDR)),], j))
      
      
      for (n in 1:nrow(tmpRes))
      {    
        system(paste("echo",tmpRes$snps[n],tmpRes$gene[n],tmpRes$statistic[n],tmpRes$pvalue[n],tmpRes$FDR[n],tmpRes$beta[n],tmpRes$myFDR[n],tmpRes$j[n], paste0(">> ",pathFinalSentinalised,tissue,"/",geneID)))
      }
      
      
      
      my.expr0 <- as.matrix( expr[[tissue]] )
      
      my.expr0 <- as.matrix( my.expr0[ ,which(colMeans( is.na(my.expr0) ) != 1) ])
      
      ## check if there is no expression for the tissue
      if(ncol(my.expr0)<1){next}
      my.markers0 <- as.matrix(markers[ , is.element(colnames(markers),rownames(my.expr0))])
      
      
      ### add the genetic covariants
      ## This is to keep the dosage covariats of the previous SNP
      if(j<2)
      { 
        my.covTMP <- read.table.rows(paste0(genotypeFile), keepRows=rownames(my.expr0), sep=" ",header=F)
        my.cov0 <- as.matrix(my.covTMP[colnames(my.markers0),2:4])
        my.cov0 <- t(my.cov0)
      }
      ## include in the covariants the dosage from the SNP choosen as "lead"
      covDose <- markers[which(SNP == rownames(markers)),colnames(my.cov0)]
      stopifnot(identical(colnames(covDose),colnames(my.cov0)))
      my.cov0 <- rbind(my.cov0,covDose)
      ##print(my.cov0)
      ## transpose needed to have the individuals as columns
      my.expr0 <- as.matrix(my.expr0[colnames(my.markers0),])
      ## in case number of columns is equal to 1 we assign again the name of the column again
      if (ncol(my.expr0)==1)
      {
        colnames(my.expr0) <- colnames(expr[[tissue]])
      }
      
      my.expr0 <- t(my.expr0)
      stopifnot( identical( colnames(my.expr0), colnames(my.markers0) ) )
      stopifnot( identical( colnames(my.cov0), colnames(my.markers0) ) )
      ## create the object SlicedData
      my.expr    <- SlicedData$new()
      my.expr$CreateFromMatrix( my.expr0 )
      my.markers <- SlicedData$new()
      my.markers$CreateFromMatrix(my.markers0)
      my.cov <- SlicedData$new()
      my.cov$CreateFromMatrix(as.matrix(my.cov0))
      rm(my.expr0, my.markers0)
      
      ## outputfile
      outputFile=paste0(tmpFolder,tissue,"/",geneID)
      ## run the tests
      store <- Matrix_eQTL_main( my.markers, my.expr, my.cov,output_file_name = NULL,pvOutputThreshold=1e-1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
      
      ## Calculate my own FDR
      
      ## we calculate the FDR without taking the min as the method used by Shabalin
      
      pval <- store$all$eqtls$pvalue
      myFDR <- sapply(pval, function(pval) p.adjust(pval, method="fdr", n =  nrow(markers.info))) 
      rm(pval)
      my.eQTLstmp <- cbind(store$all$eqtls, myFDR)
      ##print(nrow(markers.info))
      ##print(head(my.eQTLstmp)) 
      
      
      eQTLs <- my.eQTLstmp[which(my.eQTLstmp$myFDR <= FDRthr),]
      rm(my.eQTLstmp) 
      ##print(head(my.eQTLs))
      
      if (nrow(eQTLs)>0)
      {
        
        write.table(eQTLs,outputFile,row.names=F)
      }
      else
      {
        break
      }
      rm(eQTLs)    
      j=j+1
    }
    
    
  }
}


defExonicRegions <- function(regionsGene)
{
  ## regionsGene: matrix with the structure of the exons of a gene
  ## we sort the exons based on the genomic start
  if(is.unsorted(regionsGene$Exon.Chr.Start..bp.,))
  {
    regionsGene <- regionsGene[order(regionsGene$Exon.Chr.Start..bp.),]
  }
  ## we check whether there are regions with length 1bp
  if(length(which(regionsGene[,1]==regionsGene[,2]))>0)
  {  
    regionsGene <- regionsGene[-which(regionsGene[,1]==regionsGene[,2]),]
  }
  ## we check that there is no overlpapping region anin case me merge the region
  j <-2
  for (i in 2:nrow(regionsGene))
  { 
    
    if(j>nrow(regionsGene)){  
      break
    }
    
    if((regionsGene[j,1] - regionsGene[j-1,2]) <1)
    {
      ## we check whether the exon is not included already in the exon
      if (regionsGene[j,2] - regionsGene[j-1,2]>0){
        regionsGene[j-1,2] <- regionsGene[j,2]
        regionsGene <- regionsGene[-j,]
      }else{
        regionsGene <- regionsGene[-j,]
      }
      
    }else{
      j <-j+1
    }
  }
  rm(j,i)
  return(regionsGene)
}

#' @description = This function returns the length of the gene taking only into account the exonic regions
#'  
getRegionsWidth <- function(geneID,exonsdef)
{
  
  regionsGene <- exonsdef[which(exonsdef$Ensembl.Gene.ID %in% geneID),]
  ## remove duplicates
  regionsGene <- unique(regionsGene)
  ## regionsGene: matrix with the structure of the exons of a gene
  regionsGene <- defExonicRegions(regionsGene)
  ## write the GC content
  sums <- regionsGene$Exon.Chr.End..bp. -regionsGene$Exon.Chr.Start..bp.
  return(c(paste(geneID),sum(sums)))
  
}


#' @description = return the regions the regions in BED format
#'  
getRegionsBED <- function(geneID,exonsdef)
{
  
  regionsGene <- exonsdef[which(exonsdef$Ensembl.Gene.ID %in% geneID),]
  ## remove duplicates
  regionsGene <- unique(regionsGene)
  ## regionsGene: matrix with the structure of the exons of a gene
  regionsGene <- defExonicRegions(regionsGene)
  ## write the GC content
  return(as.data.frame(paste(regionsGene$Chromosome.Name,regionsGene$Exon.Chr.Start..bp.,regionsGene$Exon.Chr.End..bp.,geneID,sep='\t')))
}

#' @description This function returns the ratio of the GC exonic GC content given all the regions
#' 
ratioGCcontent <- function(geneID,GCcontentTab)
{
  ## function that calculate GC content for the exonic regions defined in the file output from GATK and it save the result in a 
  ## output file
  ## geneID: geneID from ensembl
  ## GCcontent: matrix with the GC content by exons
  
  GC <- GCcontentTab[which(GCcontentTab[,4] %in% geneID),6] 
  return(c(as.character(geneID), (sum(GC)/length(GC))))
}





