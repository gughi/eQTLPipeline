defExonicRegions <-
function(regionsGene)
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
