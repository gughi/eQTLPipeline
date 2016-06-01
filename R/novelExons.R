## These functions allow to identify new exon junctions connected between intergenic regions and nearest gene.

getOverlap <- function(query,juncRanges)
{  
  oveTmp <- findOverlaps(query,juncRanges,ignore.strand=FALSE)
  
  if(length(oveTmp) > 0)
  {
    res <- paste(as.data.frame(oveTmp)[,2],collapse = ",")
  }else{
    res <- NA
  }
  return(res)
}


getCountsJunc <- function(tmp,i,juncRanges,GR)
{
  if(!is.na(tmp$juncOverlapIndex[i]))
  {
    tmpJunc <- elementMetadata(juncRanges[as.numeric(unlist(strsplit(x = as.character(tmp$juncOverlapIndex[i]),","))),])[,2]
    overlapJunc <- overlapsAny(GR[as.character(tmp$gene[i]),],juncRanges[setdiff(which(elementMetadata(juncRanges)[,2] %in% tmpJunc),
                                                                                 (as.numeric(unlist(strsplit(x = as.character(tmp$juncOverlapIndex[i]),",")))))])
    
    counts <- sum(elementMetadata(juncRanges[setdiff(which(elementMetadata(juncRanges)[,2] %in% tmpJunc),
                                                     (as.numeric(unlist(strsplit(x = as.character(tmp$juncOverlapIndex[i]),",")))))])[,1])
    
    res <- c(as.character(overlapJunc),as.character(counts))
    rm(counts,overlapJunc,tmpJunc)
  }else{
    
    res <- c("FALSE","0")
  }
  
  return(res)
}

## This function identify whether any intergenic region has a en exon exon junction that overlaps teh intergenic region and the nearest gene 
## givinfg us infromation whther this is a novel exon
## Parameters: a table with the intergenic regions, the path of the file .bed with the junctions coming from tophat, the genomic range object with 
## definition of the genes

identifyNewExon <- function(interTable,juncPath,GR)
{
  ## load the junctions information from tophat
  junc <- read.delim(paste0(juncPath),skip = 1,header = F)
  ## give names to the columns
  colnames(junc) <- c("chr","start","end","juncName","counts","strand", "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
  ## calculate the genomics position of the exons
  ## first exon
  endFirsExon <- as.numeric(unlist(lapply(strsplit(as.character(junc$blockSizes),","),function(x){x[1]})))
  junc$endFirsExon <- junc$start+endFirsExon-1
  ## second exon
  startSecExon <- as.numeric(unlist(lapply(strsplit(as.character(junc$blockStarts),","),function(x){x[2]})))
  junc$startSecExon <- junc$start+startSecExon
  rm(startSecExon,endFirsExon)
  
  library(GenomicRanges)
  
  tab <- as.data.frame(rbind(cbind(chr=as.character(junc$chr),start=as.numeric(junc$start),end=as.numeric(junc$endFirsExon),counts=junc$counts),
                             cbind(chr=as.character(junc$chr),start=as.numeric(junc$startSecExon),end=as.numeric(junc$end),counts=junc$counts)))
  
  
  tab$name <- c(as.character(junc$juncName),as.character(junc$juncName))
  tab$start <- as.numeric(as.character(tab$start))
  tab$end <- as.numeric(as.character(tab$end))
  tab$counts <- as.numeric(as.character(tab$counts))
  
  juncRanges <- makeGRangesFromDataFrame(tab,
                                         keep.extra.columns = T )
  
  juncRanges <- sortSeqlevels(juncRanges)
  juncRanges <- sort(juncRanges)
  save(juncRanges,file=paste0(unlist(strsplit(x = as.character(juncPath),"\\."))[1],"Ranges.rda"))
  rm(tab)
  
  
  library(doParallel)
  library(foreach)
  
  detectCores()
  ## [1] 24
  # create the cluster with the functions needed to run
  cl <- makeCluster(6)
  clusterExport(cl, c("findOverlaps","paste","getOverlap","GRanges","IRanges","Rle"))
  registerDoParallel(cl)
  getDoParWorkers()
  start <- Sys.time()
  juncOverlap <- foreach(i=1:nrow(interTable),.combine=c,.verbose=F)%dopar%getOverlap(GRanges(seqnames = Rle(gsub("chr","",interTable$chr[i])),
                                                                                              ranges = IRanges(start=interTable$start[i],
                                                                                                               end = interTable$end[i])),juncRanges)
  end <- Sys.time()
  end-start
  stopCluster(cl)
  
  rm(cl,end,start)
  
  ## merge the junctions that overlap with the original table
  interTable <- cbind(interTable,juncOverlapIndex=juncOverlap)
  ## the juncOverlapIndex column tell us which are the indexes of the table juncRanges that overlap with the region if NA, no junction overlap
  head(interTable)
  
  
  ## this functions tells whether there is a junction overlapping the intergenic region and if yes, returns the counts
  rm(juncOverlap)
  
  detectCores()
  ## [1] 24
  # create the cluster with the functions needed to run
  cl <- makeCluster(6)
  clusterExport(cl, c("getCountsJunc","elementMetadata","overlapsAny","setdiff","unlist","strsplit","as.character","as.numeric","sum","which"))
  registerDoParallel(cl)
  getDoParWorkers()
  start <- Sys.time()
  JuncCounts <- foreach(i=1:nrow(interTable),.combine=rbind,.verbose=F)%dopar%getCountsJunc(interTable,i,juncRanges,GR)  
  end <- Sys.time()
  end-start
  stopCluster(cl)
  
  rm(cl,end,start,juncRanges)
  
  colnames(JuncCounts) <- c("novelExon","counts")
  
  interTable <- cbind(interTable,JuncCounts)
  
  return(interTable)
}

overlappingBP <- function(GR1,GR2)
{
  tmp <- intersect(GR1,GR2,ignore.strand=TRUE)
  if(length(tmp)>0){
    res <- abs((width(tmp) - width(GR1)))
    rm(tmp) 
  }else{
    
    res <- width(GR1)
    rm(tmp)
  }
  return(res)
}

