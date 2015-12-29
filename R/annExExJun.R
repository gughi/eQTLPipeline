## this function annotate the exon exon junctions
annExExJun <- function(exexJunID,mapExon)
{
  exExJun <- strsplit(as.character(exexJunID),"_")
  geneID <- as.character(mapExon[mapExon$exonID %in% as.character(exExJun[[1]][1]),"geneID"])
  tmp1 <- unlist(strsplit(geneID,"_"))[1]
  geneID <- as.character(mapExon[mapExon$exonID %in% as.character(exExJun[[1]][2]),"geneID"])
  tmp2 <- unlist(strsplit(geneID,"_"))[1]
  transID <- as.character(mapExon[mapExon$exonID %in% as.character(exExJun[[1]][1]),"transID"])
  transtmp1 <- unlist(strsplit(transID,"_"))[2]
  transID <- as.character(mapExon[mapExon$exonID %in% as.character(exExJun[[1]][2]),"transID"])
  transtmp2 <- unlist(strsplit(transID,"_"))[2]

  ##  stopifnot(identical(tmp1,tmp2))
  if(!identical(tmp1,tmp2))
  {
    tmp1 <- paste(tmp1,tmp2,sep=";")  
  }
  
  if(!identical(transtmp1,transtmp2))
  {
    transtmp1 <- paste(transtmp1,transtmp2,sep=";")  
  }
  
  exon1 <- mapExon[which(mapExon$exonID %in% as.character(exExJun[[1]][1])),c("chr","start","end")]
  exon2 <- mapExon[which(mapExon$exonID %in% as.character(exExJun[[1]][2])),c("chr","start","end")]
  coor <- c(tmp1,transtmp1,exon1,exon2)
  rm(tmp1,tmp2,exon1,exon2)
  names(coor) <- c("geneID","transID", "chrExon1","startExon1","endExon1","chrExon2","startExon2","endExon2")
  coor <- as.data.frame(coor)
  rownames(coor) <- as.character(exexJunID)
  return(coor)  

}

## the function below returns the number of the ex-ex junctions from the same Gene targeted by the same SNP 
# 
numOfSameTranscript <- function(SNP,gene,eQTLlist)
{
  tmp<-eQTLlist[which(eQTLlist$snps %in% SNP),]
  return(table(tmp$geneID)[gene])  
  
}


## this function check whether the exon does not overlap with any other region
uniqueExon <- function(chr,exon1Start,exon1End,exon2Start,exon2End,ensembl)
{
  
  library(biomaRt)
  library(Gviz)
  
  ##geneType <- "lincRNA"
  
  ##gene <- "ENSG00000006555"
  
  chromoReg <- paste0(chr,":",exon1Start,":",exon2End,":-1,",
                      chr,":",exon1Start,":",exon2End,":1")
  filters <- c("chromosomal_region")
  defGen <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                               "exon_chrom_end","strand","gene_biotype",
                               "ensembl_gene_id",'ensembl_exon_id',
                               'ensembl_transcript_id',"external_gene_id"),
                  filters=filters, values=list(chromoReg), mart=ensembl)
  
  GR <- GRanges(seqnames = Rle(chr),
                ranges = IRanges(start=c(exon1Start,exon2Start),
                                 end = c(exon1End,exon2End),
                                 names = c("exon1","exon2")))
  
  counts<- countOverlaps(GR, GRanges(seqnames = Rle(defGen$chromosome_name),
                                   ranges = IRanges(start=defGen$exon_chrom_start,
                                                    end = defGen$exon_chrom_end,
                                                    names = defGen$ensembl_exon_id)))
  # removed the exon it self
  counts<-counts-1
  if(sum(counts)<2){
     return(TRUE)
  }else if(counts[1]==0 || counts[2]==0){
    return(TRUE)
  }else{
    return(FALSE)
  }
    
  
}

## return the exonrank for the transcript give a chr, start and stop position
## used for annotate the exon exon junctions
getRank <- function(chr,start,end,rankTab)
{  
  GR <- GRanges(seqnames = Rle(rankTab$chromosome_name),
                ranges = IRanges(start=rankTab$exon_chrom_start,
                                 end = rankTab$exon_chrom_end),
                rank=rankTab$rank)
  
  exGR <- GRanges(seqnames = Rle(chr),
                  ranges = IRanges(start=start,
                                   end = end))
  
  res <- subsetByOverlaps(GR,exGR)
  return(paste(res$rank,sep=","))
}



