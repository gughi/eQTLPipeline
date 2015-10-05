# exons is a vector of two values containing the IDs of the two exons that form the juction
GCcontJunction <- function(exons,GCcontJunc)
{
  GC1 <- GCcontentTab[which(GCcontentTab$X4_usercol %in% exons[1]), "X6_pct_gc"]   
  GC2 <- GCcontentTab[which(GCcontentTab$X4_usercol %in% exons[2]), "X6_pct_gc"]
  return((GC1+GC2)/2)
}
  