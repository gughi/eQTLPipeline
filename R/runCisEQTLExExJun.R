runCisEQTLGene <- function(i,resids,snpLocation,outputFolder,my.covTMP,fullResults,MySQL=FALSE)
{
  
  exExJunID  <- colnames(resids)[i]
  expr <-  as.matrix(resids[,i]); rm(i)
  colnames(expr) <- as.character(exExJunID)
  print(exExJunID)
  
  ## load the snp data
  geneID <- as.character(mapExon[mapExon$exonID %in% unlist(strsplit(exExJunID,"_"))[1],"geneID"])
  geneID <- unlist(strsplit(geneID,"_"))[1]
  fn.rda <- paste0(snpLocation,"/", geneID, ".rda")
  load(fn.rda)
  rm(fn.rda)
  
  ## create folder for the output
  dir.create(file.path(outputFolder),showWarnings=FALSE)
  dir.create(file.path(fullResults),showWarnings=FALSE)
  
  
  if(markers=="No polymorphisms found within +/- 1Mb of TSS or TES")
  {
    stop("No polymorphisms found within +/- 1Mb of TSS or TES")
  }
  
  library(MatrixEQTL)
  
  
  ## check if expr and markers have the same number of individuals, basically if all the individuals 
  ## have expression data for that gene but also markers
  
  expr <- as.matrix(expr[ , which(colMeans( is.na(expr) ) != 1) ])
  
  my.markers0 <- as.matrix(markers[ , is.element(colnames(markers),rownames(expr))])
  
  ### add the genetic covariants
  
  my.cov0 <- as.matrix(my.covTMP[colnames(my.markers0),2:4])
  my.cov0 <- t(my.cov0)
  
  
  ## transpose needed to have the individuals as columns
  my.expr0 <- as.matrix(expr[colnames(my.markers0),])
  ## in case number of columns is equal to 1 we assign again the name of the column again
  if (ncol(my.expr0)==1)
  {
    colnames(my.expr0) <- geneID
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
  outputFile=paste0(outputFolder,geneID)
  store <- Matrix_eQTL_main( my.markers, my.expr, my.cov,output_file_name = paste0(fullResults,geneID),pvOutputThreshold=1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
  
  ## we calculate the FDR without taking the min as the method used by Shabalin
  pval <- store$all$eqtls$pvalue
  myFDR <- sapply(pval, function(pval) p.adjust(pval, method="fdr", n =  nrow(markers.info))) 
  rm(pval)
  my.eQTLstmp <- cbind(store$all$eqtls, myFDR)
  
  ## we save it in the database
  if(MySQL==TRUE)
  {
    
    #       con <- dbConnect(MySQL(), user="xxx", password="xxx", dbname="test", host="localhost", client.flag=CLIENT_MULTI_STATEMENTS)
    #       dbWriteTable(conn=con,name=as.character(tissue),value=my.eQTLstmp,row.names=F,append=T)
    #       dbDisconnect(con)
  }
  
  my.eQTLs <- my.eQTLstmp[which(my.eQTLstmp$myFDR <= 0.10),]
  rm(my.eQTLstmp) 
  ##print(head(my.eQTLs))
  
  print(nrow(my.eQTLs))
  if (nrow(my.eQTLs)>0)
  {
    write.table(my.eQTLs,outputFile,row.names=F)
  }
  
  rm(my.expr,my.markers,my.cov, markers.info,my.eQTLs,store,myFDR)
}

