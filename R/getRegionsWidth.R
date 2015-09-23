getRegionsWidth <-
function(geneID,exonsdef,intronic=FALSE)
{
  regionsGene <- exonsdef[which(exonsdef$Ensembl.Gene.ID %in% geneID),]
  ## remove duplicates
  regionsGene <- unique(regionsGene)
  ## regionsGene: matrix with the structure of the exons of a gene
  if(!intronic)
  {
    regionsGene <- defExonicRegions(regionsGene)
  }
  ## write the GC content
  sums <- regionsGene$Exon.Chr.End..bp. -regionsGene$Exon.Chr.Start..bp.
  return(c(paste(geneID),sum(sums)))
}
