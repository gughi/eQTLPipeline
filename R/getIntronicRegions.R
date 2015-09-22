

getIntronicRegions <- function(geneID,exonsdef)
{
  ## geneID:: the gene ensembl ID we want to define the intronic regions
  ## exonsdef: the matrix that defines the exonic regions
  ## it return a BED format of the intonic regions defined
  regionsGene <- exonsdef[which(exonsdef$Ensembl.Gene.ID %in% geneID),]
  ## remove duplicates
  regionsGene <- unique(regionsGene)
  ## regionsGene: matrix with the structure of the exons of a gene
  regionsGene <- defExonicRegions(regionsGene)
  
  intronicRegionsGene <- NULL
  if(nrow(regionsGene)==1)
  {
    intronicRegionsGene <- t(as.data.frame(c("NA","NA",as.character(regionsGene[,4]),"NA")))   
  }else{
    for(i in 1:(nrow(regionsGene)-1))
    {
      intronicRegionsGene <- rbind(intronicRegionsGene,c(regionsGene[i,2],regionsGene[i+1,1],as.character(regionsGene[i,4]),as.character(regionsGene[i,5])))
    }  
    
    if(length(which(intronicRegionsGene[,1]==intronicRegionsGene[,2]))>0)
    {  
      intronicRegionsGene <- intronicRegionsGene[-which(intronicRegionsGene[,1]==intronicRegionsGene[,2]),]
    }
    
  }
  ## we check whether there are regions with length 1bp
  
  return(as.data.frame(paste(intronicRegionsGene[,4],intronicRegionsGene[,1],intronicRegionsGene[,2],geneID,sep='\t')))
  
}

