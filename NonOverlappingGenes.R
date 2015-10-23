## script top find genomic overlaps between genes

library(GenomicRanges)

genesMap <- read.delim("data/general/ensemblGenes.txt")

head(genesMap)

## define the genomic regions
GR <- GRanges(seqnames = Rle(genesMap$Chromosome.Name),
              ranges = IRanges(start=genesMap$Gene.Start..bp.,
                               end = genesMap$Gene.End..bp.,
                               names = genesMap$Ensembl.Gene.ID))





library(doParallel)
library(foreach)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(20)
clusterExport(cl, c("overlapsAny"))

registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
nonOveGenes <- foreach(i=1:length(GR[2,]),.combine=c,.verbose=F)%dopar%overlapsAny(GR[i,],GR[-i,],ignore.strand=TRUE)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)








