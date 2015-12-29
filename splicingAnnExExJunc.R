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


# exExJunAnn <- exExJunAnn[which(exExJunAnn$geneID %in% as.character(geneNames$ensembl_gene_id)),]

## we get the table for the exons IDs
#load("data/results/finaleQTLs/exExJun.unsentinalised.rda")
load("data/expr/rawCounts/genic/fullExExJun.rda")
rm(map)

library(doParallel)
library(foreach)
library(GenomicRanges)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(5)
clusterExport(cl, c("grep"))
registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
exonIDs <- foreach(i=1:length(unique(geneNames$ensembl_gene_id)),.verbose=F)%dopar%(mapExon[grep(pattern=as.character(unique(geneNames$ensembl_gene_id)[i]),
                                                                                                                      x=as.character(mapExon$geneID)),"exonID"])
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)
rm(end,start)
exonIDs <- unlist(exonIDs)

## get the exon-Exon junctions expressed
idx <- which(as.character(expr$Exon1ID) %in% as.character(exonIDs))
idxtmp <- which(as.character(expr$Exon2ID) %in% as.character(exonIDs))
idx <- unique(c(idx,idxtmp))
rm(idxtmp,exonIDs)
exExIDs <- paste(expr[idx,"Exon1ID"],expr[idx,"Exon2ID"],sep="_")
rm(idx)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(5)
clusterExport(cl, c("annExExJun"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
exExJunAnn <- foreach(i=1:length(exExIDs),.combine=rbind,.verbose=F)%dopar%(annExExJun(exExIDs[i],mapExon))
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)
rm(end,start,mapExon)

save(exExJunAnn,file="data/general/neuroGenesAnnExExJun.PUTM.rda")


detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(10)
clusterExport(cl, c("getRank","Rle","IRanges","subsetByOverlaps","GRanges"))
registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
rankExonOne <- foreach(i=1:nrow(exExJunAnn),.combine=c,.verbose=F)%dopar%getRank(as.character(exExJunAnn[i,"chrExon1"]),
                                                        as.numeric(as.character(exExJunAnn[i,"startExon1"])),
                                                        as.numeric(as.character(exExJunAnn[i,"endExon1"])),
                                                        geneNames[which(as.character(geneNames$ensembl_transcript_id) %in% as.character(exExJunAnn[i,"transID"])),])

rankExonTwo <- foreach(i=1:nrow(exExJunAnn),.combine=c,.verbose=F)%dopar%getRank(as.character(exExJunAnn[i,"chrExon2"]),
                                                                                 as.numeric(as.character(exExJunAnn[i,"startExon2"])),
                                                                                 as.numeric(as.character(exExJunAnn[i,"endExon2"])),
                                                                                 geneNames[which(as.character(geneNames$ensembl_transcript_id) %in% as.character(exExJunAnn[i,"transID"])),])


##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)

exExJunAnn$rankFirstExon <- rankExonOne
exExJunAnn$rankSecondExon <- rankExonTwo


load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.PUTM")

expr.RPKM.cqn <- RPKM.cqn[intersect(rownames(RPKM.cqn),rownames(exExJunAnn)),]

save(exExJunAnn,expr.RPKM.cqn,file="data/general/neuroGenesExprAndAnn.PUTM.rda")

head(exExJunAnn,1000)
