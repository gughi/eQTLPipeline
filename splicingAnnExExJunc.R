##Annotation Novel exon-exon Junctions

gene <- read.delim("data/general/NP4.6c.raw.n.bed")
gene <- unique(gene$LOCUS)

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id",
                                "chromosome_name","ensembl_transcript_id",
                                "gene_biotype","exon_chrom_start","exon_chrom_end","rank"),
                   verbose = T,
                   filters="hgnc_symbol",
                   values=gene, mart=ensembl)

geneNames <- geneNames[-which(as.character(geneNames$gene_biotype) %in% "LRG_gene"),]
load("data/general/annotationExExJun.PUTM.rda")
rm(gene,ensembl)


exExJunAnn <- exExJunAnn[which(exExJunAnn$geneID %in% as.character(geneNames$ensembl_gene_id)),]

## we get the table for the exons IDs
head(exExJunAnn)


head(exExJunAnn)


head(geneNames)


singleExExJun <- exExJunAnn[1,]
chr <- singleExExJun$chrExon1
start <- singleExExJun$startExon1
end <- singleExExJun$endExon1
head(singleExExJun)

rankTab <- geneNames[which(as.character(geneNames$ensembl_transcript_id) %in% as.character(singleExExJun$transID)),]

head(rankTab)


getRank(chr,start,end, rankTab)



library(doParallel)
library(foreach)
library(GenomicRanges)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(10)
clusterExport(cl, c("getRank","Rle","IRanges","subsetByOverlaps","GRanges"))
registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
rankExonOne <- foreach(i=1:10,.verbose=F)%dopar%getRank(exExJunAnn[i,"chrExon1"],
                                                        exExJunAnn[i,"startExon1"],
                                                        exExJunAnn[i,"endExon1"],
  geneNames[which(as.character(geneNames$ensembl_transcript_id) %in% as.character(exExJunAnn[i,"transID"])),])
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)

dim(unique(exExJunAnn))




