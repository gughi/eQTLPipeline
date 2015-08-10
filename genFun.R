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
#' @param exprThe expression matrix
#' @param covariates to use for the correction
#' @param covs the covariates to exclude from PEER analysis
#' @param outputFile where the results are going to be saved
#' 
doResidualCorrection <- function(expr,covs,outputFile)
{
  ## outputFolder <- "/home/seb/expressionData/exonExpr/"
  ## covs is the PEER
  ## load expression data un filtered
  ## expr is matrix nxm where n are the samples and m the genes
  print("Loading the expression")
  
  ## load the knowing factors (age and sex ) and unknown factors (PEER axes)
  ## Residual correction with linear model, step 3
  print("Performing the residual correction")
  resids <- apply(exprAll, 2, function(y){
    lm( y ~ . , data=covs)$residuals 
  })
  
  print("Saving the residuals")
  save(resids,paste0(outputFolder,"/resCorExpr13.csv"))
  resids
}

