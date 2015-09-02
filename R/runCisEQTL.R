runCisEQTL <-
function(i,ensemblGenes,exprLocation,snpLocation,outputFolder,genotypeFile,fullResults,MySQL=FALSE)
{
  library(RMySQL) 
  
  
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

    ## we save it in the database
    if(MySQL==TRUE)
    {
      con <- dbConnect(MySQL(), user="adai", password="adai123", dbname="test", host="localhost", client.flag=CLIENT_MULTI_STATEMENTS)
      dbWriteTable(conn=con,name=as.character(tissue),value=my.eQTLstmp,row.names=F,append=T)
      dbDisconnect(con)
    }
    else
    {
      save(my.eQTLstmp,file=paste0(fullResults,tissue,"/",geneID,".rda"),compress="bzip2")
    }
    
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
