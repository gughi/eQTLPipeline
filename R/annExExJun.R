## this function annotate the exon exon junctions
annExExJun <- function(exexJunID,mapExon)
{
  exExJun <- strsplit(as.character(exexJunID),"_")
  geneID <- as.character(mapExon[mapExon$exonID %in% as.character(exExJun[[1]][1]),"geneID"])
  tmp1 <- unlist(strsplit(geneID,"_"))[1]
  geneID <- as.character(mapExon[mapExon$exonID %in% as.character(exExJun[[1]][2]),"geneID"])
  tmp2 <- unlist(strsplit(geneID,"_"))[1]
##  stopifnot(identical(tmp1,tmp2))
  if(!identical(tmp1,tmp2))
  {
    tmp1 <- paste(tmp1,tmp2,sep=";")  
  }
  exon1 <- mapExon[which(mapExon$exonID %in% as.character(exExJun[[1]][1])),c("chr","start","end")]
  exon2 <- mapExon[which(mapExon$exonID %in% as.character(exExJun[[1]][2])),c("chr","start","end")]
  coor <- c(tmp1,exon1,exon2)
  rm(tmp1,tmp2,exon1,exon2)
  names(coor) <- c("geneID", "chrExon1","startExon1","endExon1","chrExon2","startExon2","endExon2")
  coor <- as.data.frame(coor)
  rownames(coor) <- as.character(exexJunID)
  return(coor)
}
