

library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

## in this scrip t we identify if there is any junction in our data for the 
load(file="data/general/defIntergenic.rda")

## we remove Intergenic regions that don't have any new portion but that are totally inside the gene

tmp <- regInformation[which(regInformation$uniquePortion > 0),]
tmp$overlap <- as.numeric(as.character(tmp$overlap))
tmp <- tmp[which(tmp$overlap < 2),]

## we load the regions of the exon junction ( it means exon and junctions inclueded) from tophat
junc <- read.delim("data/junctions/junctions.bed",skip = 1,header = F)

colnames(junc) <- c("chr","start","end","juncName","counts","strand", "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")

endFirsExon <- as.numeric(unlist(lapply(strsplit(as.character(junc$blockSizes),","),function(x){x[1]})))
junc$endFirsExon <- junc$start+endFirsExon-1

startSecExon <- as.numeric(unlist(lapply(strsplit(as.character(junc$blockStarts),","),function(x){x[2]})))
junc$startSecExon <- junc$start+startSecExon

rm(startSecExon,endFirsExon)

library(GenomicRanges)

tab <- as.data.frame(rbind(cbind(chr=as.character(junc$chr),start=as.numeric(junc$start),end=as.numeric(junc$endFirsExon),counts=junc$counts),
                                cbind(chr=as.character(junc$chr),start=as.numeric(junc$startSecExon),end=as.numeric(junc$end),counts=junc$counts)))

tail(tab)

tab$name <- c(as.character(junc$juncName),as.character(junc$juncName))
tab$start <- as.numeric(as.character(tab$start))
tab$end <- as.numeric(as.character(tab$end))
tab$counts <- as.numeric(as.character(tab$counts))

juncRanges <- makeGRangesFromDataFrame(tab,
                         keep.extra.columns = T )


juncRanges <- sortSeqlevels(juncRanges)
juncRanges <- sort(juncRanges)

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


tmp <- cbind(tmp,juncOverlapIndex=juncOverlap)


# ## examples
# tmp["DER362",]
# juncRanges[3762,]
# 
# ## example of an extention of the exon
# tmp["DER276",]
# juncRanges[2773,]
# 
# 
# tmp["DER592",]
# juncRanges[9031,]
# 
# 
# 
# tmp["DER1395",]
# juncRanges[19852,]



load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

## 5%FDR
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.05),1:7]
head(eQTLPUTM)



length(intersect(rownames(tmp[which(!is.na(juncOverlap)),]),eQTLPUTM$gene))

## we look at one example with eQTLs

head(intersect(rownames(tmp[which(!is.na(juncOverlap)),]),eQTLPUTM$gene),20)

head(tmp[intersect(rownames(tmp[which(!is.na(juncOverlap)),]),eQTLPUTM$gene),],30)


## now we will check whether the intergenic region identified is an new exon


eQTLsIntergenic <- tmp[intersect(rownames(tmp[which(!is.na(juncOverlap)),]),eQTLPUTM$gene),]

## we select the regions that have distance at from the gene at least > 0bp
eQTLsIntergenic <- eQTLsIntergenic[which(as.numeric(as.character(eQTLsIntergenic$distance)) >0),]


head(eQTLsIntergenic)

genesMap <- read.delim("data/general/ensemblGenes.txt")

head(genesMap)

## define the genomic regions
GR <- GRanges(seqnames = Rle(genesMap$Chromosome.Name),
              ranges = IRanges(start=genesMap$Gene.Start..bp.,
                               end = genesMap$Gene.End..bp.,
                               names = genesMap$Ensembl.Gene.ID))




for(i in 1:nrow(eQTLsIntergenic))
{
  
  
    
  tmpJunc <- elementMetadata(juncRanges[as.numeric(unlist(strsplit(x = as.character(eQTLsIntergenic$juncOverlapIndex[i]),","))),])[,2]
  print(overlapsAny(GR[as.character(eQTLsIntergenic$gene[i]),],juncRanges[setdiff(which(elementMetadata(juncRanges)[,2] %in% tmpJunc),
                                       (as.numeric(unlist(strsplit(x = as.character(eQTLsIntergenic$juncOverlapIndex[i]),",")))))]))

}
























