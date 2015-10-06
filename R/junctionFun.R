# exons is a vector of two values containing the IDs of the two exons that form the juction
GCcontJunction <- function(exons,GCcontJunc)
{
  GC1 <- GCcontentTab[which(GCcontentTab$X4_usercol %in% exons[1]), "X6_pct_gc"]   
  GC2 <- GCcontentTab[which(GCcontentTab$X4_usercol %in% exons[2]), "X6_pct_gc"]
  return((GC1+GC2)/2)
}


lengthJunction <- function(exons,mapExon)
{
  Le1 <- mapExon[which(mapExon$exonID %in% exons[1]), "length"]   
  Le2 <- mapExon[which(mapExon$exonID %in% exons[2]), "length"]
  return(Le1+Le2)
}
