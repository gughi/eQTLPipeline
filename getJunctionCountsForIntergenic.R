## scripts that checks whether the intergenic regions identified has any exon-exon junction connected to the nearest gene.
## it returns a matrix with boolean value that tells whether it has or not the exon-exon junction and the counts



library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

## functions needed to run the script/function
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



## load teh intergenic regions
load(file="data/general/defIntergenic.rda")

## we remove Intergenic regions that don't have any new portion but that are totally inside the gene

tmp <- regInformation[which(regInformation$uniquePortion > 0),]
tmp$overlap <- as.numeric(as.character(tmp$overlap))
tmp <- tmp[which(tmp$overlap < 2),]



  ## load the junctions information from tophat
  junc <- read.delim("data/junctions/junctions.bed",skip = 1,header = F)
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
  juncOverlap <- foreach(i=1:nrow(tmp),.combine=c,.verbose=F)%dopar%getOverlap(GRanges(seqnames = Rle(gsub("chr","",tmp$chr[i])),
                                                                                       ranges = IRanges(start=tmp$start[i],
                                                                                                        end = tmp$end[i])),juncRanges)
  end <- Sys.time()
  end-start
  stopCluster(cl)
  
  rm(cl,end,start)
  
  ## merge the junctions that overlap with the original table
  tmp <- cbind(tmp,juncOverlapIndex=juncOverlap)
  ## the juncOverlapIndex column tell us which are the indexes of the table juncRanges that overlap with the region if NA, no junction overlap
  head(tmp)
  
  ## we select the regions that have distance at from the gene at least > 0bp
   <- eQTLsIntergenic[which(as.numeric(as.character(eQTLsIntergenic$distance)) >0),]
  
  
  genesMap <- read.delim("data/general/ensemblGenes.txt")
  
  head(genesMap)
  
  ## define the genomic regions
  GR <- GRanges(seqnames = Rle(genesMap$Chromosome.Name),
                ranges = IRanges(start=genesMap$Gene.Start..bp.,
                                 end = genesMap$Gene.End..bp.,
                                 names = genesMap$Ensembl.Gene.ID))
  

  ## this functions tells whether there is a junction overlapping the intergenic region and if yes, returns the counts
  getCountsJunc <- function(tmp,i,juncRanges)
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



  getCountsJunc(tmp,i,juncRanges)  
  
  


