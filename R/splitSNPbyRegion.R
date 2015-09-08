### Manuel Sebastian Guelfi
### 30-06-15
### This scripts put all the SNPs in cis with all the regions

splitSNPsByRegion <- function(i,allMarkers,intergenicRegions,outputFolder,logFolder,snps.map,imputed.info,regIDsLogNo){
  
  ## we select the chromosome
  ##chr <- "21
  # sys.source("/home/adai/scripts/common_functions.R",
  #             attach(NULL, name="myenv"))
  ## where the snps are going to be saved
  fn.rda <- paste0(outputFolder, rownames(intergenicRegions)[i], ".rda")
  cmd <- paste0("awk '$1==\"", gsub("chr","",intergenicRegions[i,]$chr), "\" && $2 >=",
                intergenicRegions[i,]$start - 1000000,
                " && $2 <=", intergenicRegions[i,]$end + 1000000,
                " {print $3}' ",snps.map)
  
  
  vars <- system(cmd, intern=T)
  
    if( length(vars) == 0 ){
    markers.info <- "No polymorphisms found within +/- 1Mb of TSS or TES"
    markers      <- "No polymorphisms found within +/- 1Mb of TSS or TES"
    
    save(markers, markers.info,
         file=fn.rda, compress="bzip2", ascii=T)
    
    ##system( paste("touch", fn.pval) )
    
    cat(rownames(intergenicRegions)[i], intergenicRegions[i,]$chr, intergenicRegions[i,]$start,intergenicRegions[i,]$end, 0, 0, "\n",
        file=paste(logFolder,rownames(intergenicRegions)[i],".log", sep="")) 
    cat(rownames(intergenicRegions)[i], "\n", file=paste0(regIDsLogNo), append=T)
    
    stop("Problem No Polymorphisms")
  }
  
  markers.info <- read.table.rows(paste0(imputed.info), keepRows=vars, sep=" ")
  
  ## add the type of the marker, e.g. SNP, indel
  markers.info$type <- ifelse( sapply(strsplit(rownames(markers.info), split=":"), length)==2, "SNP", "indel" )
  
  ## add to the rownames in the markers.info the suffix 'chr'
  rownames(markers.info) <- paste0("chr", rownames(markers.info))
  
  ### head(markers.info)
  ### OUTPUT:
  ##              Al1 Al2   Freq1     MAF AvgCall     Rsq        rsid genotyped type
  ##chr21:17905657   T   C 0.94402 0.05598 0.99999 1.00360  rs17276117         1  SNP
  ##chr21:17910459   A   C 0.92164 0.07836 0.99999 1.00365  rs17276124         1  SNP
  ##chr21:17912496   T   C 0.94402 0.05598 0.99999 1.00349  rs17276131         1  SNP
  ##chr21:17917427   G   A 0.92177 0.07823 0.99987 1.00202  rs17210183         1  SNP
  ##chr21:17929155   C   T 0.89289 0.10711 0.96604 0.67502 rs141445732         0  SNP
  ##chr21:17930405   C   T 0.88505 0.11495 0.99875 0.99191  rs11702768         1  SNP
  ##
  
  ### select the file to get genotyped dosages for SNPs and indels
  
  
  ## Select the marker according with vars
  markers <- allMarkers[vars,]
  
  ## head(markers)
  ##             Al1 Al2   Freq1     MAF AvgCall     Rsq 00_38 01_37 02_06 03_28 04_05 04_19 05_10 06_02 06_53 07_28 07_37 07_73
  ## 21:17905657   T   C 0.94402 0.05598 0.99999 1.00360 2.000 2.000 2.000 1.000 2.000 2.000 2.000 2.000 2.000 2.000 1.000 2.000
  ## 21:17910459   A   C 0.92164 0.07836 0.99999 1.00365 2.000 2.000 2.000 1.000 1.000 2.000 2.000 2.000 2.000 2.000 1.000 2.000
  ## 21:17912496   T   C 0.94402 0.05598 0.99999 1.00349 2.000 2.000 2.000 1.000 1.999 2.000 2.000 2.000 2.000 2.000 1.000 2.000
  ## 21:17917427   G   A 0.92177 0.07823 0.99987 1.00202 2.000 2.000 2.000 1.000 1.000 2.000 2.000 2.000 2.000 2.000 1.000 2.000
  ## 21:17929155   C   T 0.89289 0.10711 0.96604 0.67502 1.952 1.951 1.976 0.994 0.982 1.645 1.972 1.928 1.934 1.972 1.045 1.914
  ## 21:17930405   C   T 0.88505 0.11495 0.99875 0.99191 1.999 1.999 1.056 2.000 2.000 1.999 1.000 2.000 2.000 2.000 1.000 2.000
  
  ## add to the rownames in the markers.info the suffix 'chr'
  rownames(markers) <- paste0("chr", rownames(markers))
  
  ## check if the markers have been selected properly otherwise it stops adn produce an message error
  stopifnot( identical( rownames(markers), rownames(markers.info) ) )
  
  save(markers, markers.info,
       file=fn.rda, compress="bzip2", ascii=T)
  rm(fn.rda)
  
  cat(rownames(intergenicRegions)[i], intergenicRegions[i,]$chr, intergenicRegions[i,]$start,intergenicRegions[i,]$end, sum(markers.info$type=="SNP"), sum(markers.info$type=="indel"), "\n",
      file=paste(logFolder,rownames(intergenicRegions)[i],".log", sep=""))
}







