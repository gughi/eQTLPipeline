getRegionsBED <-
function(geneID,exonsdef)
{
  
  regionsGene <- exonsdef[which(exonsdef$Ensembl.Gene.ID %in% geneID),]
  ## remove duplicates
  regionsGene <- unique(regionsGene)
  ## regionsGene: matrix with the structure of the exons of a gene
  regionsGene <- defExonicRegions(regionsGene)
  ## write the GC content
  return(as.data.frame(paste(regionsGene$Chromosome.Name,regionsGene$Exon.Chr.Start..bp.,regionsGene$Exon.Chr.End..bp.,geneID,sep='\t')))
}
