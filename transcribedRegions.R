## in this script we will get the 
load("data/expr/rawCounts/genic/fullExExJun.rda")
gene <- read.delim("data/general/NP4.6c.raw.n.bed")

gene <- unique(gene$LOCUS)

library(biomaRt)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

head(listAttributes(ensembl))
grep("name",head(listFilters(ensembl),300))

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","chromosome_name","start_position","end_position","gene_biotype"),
                   verbose = T,
                   filters="hgnc_symbol",
                   values=gene, mart=ensembl)

## remove all the LRG genes
geneNames <- geneNames[-which(geneNames$gene_biotype %in% "LRG_gene"),]

geneNames <- geneNames[order(geneNames$chromosome_name),]


load("data/general/overlappingGenes.rda")

neuroNonOveGen <- geneNames[-which(geneNames[,1] %in% rownames(as.data.frame(listNonOve))),]

library("devtools")
load_all()

library(doParallel)
library(foreach)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(5)
clusterExport(cl, c("novelTransRegion","getBM","subsetByOverlaps"))
registerDoParallel(cl)
getDoParWorkers()

## remove gene that comes from the patch HG987_PATCH - gene KCNJ12
neuroNonOveGen <- neuroNonOveGen[c(-58-),]

neuroNonOveGen <- neuroNonOveGen[c(-59,-61,-62,63),]

start <- Sys.time()
novelRegions <- foreach(i=1:nrow(neuroNonOveGen),.combine=rbind,.verbose=F)%dopar%novelTransRegion(neuroNonOveGen[i,],ensembl)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)

sum(unlist(novelRegions[,3]),na.rm=T)



load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
IDs=PUTM$A.CEL_file

plotReadDepth(gene=neuroNonOveGen[i,1],ensembl=ensembl,IDs=IDs)

ensemblv75 <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Feb2014.archive.ensembl.org",
                      dataset="hsapiens_gene_ensembl")
exonDef <- getBM(attributes=c("ensembl_gene_id","chromosome_name","ensembl_exon_id","exon_chrom_start","exon_chrom_end"),
                 verbose = T,
                 filters="ensembl_gene_id",
                 values=neuroNonOveGen[i,1], mart=ensemblv75)

exonDef <- GRanges(paste0("chr",exonDef[,2]), IRanges(exonDef[,4], exonDef[,5]))
## select the transcribed regions identified in the data
tmp <- subsetByOverlaps(expressedRegions$chr1$regions,
                        GRanges(paste0("chr",neuroNonOveGen[i,3]),
                                IRanges(neuroNonOveGen[i,4], neuroNonOveGen[i,5])))
nrow(as.data.frame(tmp))

table(countOverlaps(tmp[which(tmp$value>10),], exonDef)==0)

plotReadDepth(gene=neuroNonOveGen[4,1],ensembl=ensembl,IDs=IDs)





