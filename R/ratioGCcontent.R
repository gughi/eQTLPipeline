#' @description This function returns the ratio of the GC exonic GC content given all the regions
ratioGCcontent <-
function(geneID,GCcontentTab)
{
  ## function that calculate GC content for the exonic regions defined in the file output from GATK and it save the result in a 
  ## output file
  ## geneID: geneID from ensembl
  ## GCcontent: matrix with the GC content by exons
  
  GC <- GCcontentTab[which(GCcontentTab[,4] %in% geneID),6] 
  return(c(as.character(geneID), (sum(GC)/length(GC))))
}


ratioGCcontentExonExon <-
  function(juncID,GCcontentTab)
  {
    ## function that calculate GC content for the  
    ## output file
    ## geneID: geneID from ensembl
    ## GCcontent: matrix with the GC content by exons
    
    exons <- unlist(strsplit(juncID,"_"))
    GC <- ((GCcontentTab[which(GCcontentTab[,4] %in% exons[1]),6] + GCcontentTab[which(GCcontentTab[,4] %in% exons[2]),6])/2)
    return(c(as.character(juncID), (sum(GC)/length(GC))))
  }


