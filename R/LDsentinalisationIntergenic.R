LDsentinalisationIntergenic <-
function(exprIntergenic,regID,pathFinalSentinalised,pathUnsentinalised,FDRthr,my.covTMP,snpLocation,tmpFolder)
{
  
  library(MatrixEQTL)
  expr <- as.matrix(exprIntergenic[,as.character(regID)])
  colnames(expr) <- regID

  
  ## initialise my.cov0
  fn.rda <- paste0(snpLocation, regID, ".rda")
  load(fn.rda)
  
  rm(fn.rda)


    eQTLs <- read.delim(paste0(pathUnsentinalised,"/",regID),sep=" ")
    ## j is the degree subsignals
    j=1
    my.cov0 <- NA
    while(TRUE)
    {
      
      ## check if expr and markers have the same number of individuals, basically if all the individuals have expression data for that gene but also markers
      ## check whether there is any other eQTL
      if (j>1)
      { 
        if(!file.exists(paste0(tmpFolder,"/",regID)))
        {
          break
        }
        eQTLs <- read.delim(paste0(tmpFolder,"/",regID),sep=" ")
        
      }
      else
      {
        system(paste("echo","snps","gene","t-stat","pvalue","FDR","beta","myFDR","degree", paste0(">> ",pathFinalSentinalised,"/",regID)))
      }
      
      if (min(eQTLs$myFDR)>FDRthr)
      {
        print(min(eQTLs$myFDR))
        break
      }
      
      ##if(eQTLs$snps=="No significant associations were found." | nrow(eQTLs)<1)
      ##{
      if(file.exists(paste0(tmpFolder,"/",regID)))
      {
        system(paste0("rm ",tmpFolder,"/",regID))
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
        system(paste("echo",tmpRes$snps[n],tmpRes$gene[n],tmpRes$statistic[n],tmpRes$pvalue[n],tmpRes$FDR[n],tmpRes$beta[n],tmpRes$myFDR[n],tmpRes$j[n], paste0(">> ",pathFinalSentinalised,"/",regID)))
      }
      
      
      
      my.expr0 <- as.matrix(expr)
      
      my.expr0 <- as.matrix( my.expr0[ ,which(colMeans( is.na(my.expr0) ) != 1) ])
      
      ## check if there is no expression for the tissue
      if(ncol(my.expr0)<1){next}
      my.markers0 <- as.matrix(markers[ , is.element(colnames(markers),rownames(my.expr0))])
      
      
      ### add the genetic covariants
      ## This is to keep the dosage covariats of the previous SNP
      if(j<2)
      { 
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
        colnames(my.expr0) <- regID
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
      outputFile=paste0(tmpFolder,"/",regID)
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
