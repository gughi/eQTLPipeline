## in this script we will get the 
## load("data/expr/rawCounts/genic/fullExExJun.rda")

## load libraries and data
library(devtools)
setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()
library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

## load neuro genes passed by Jana
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
## load the genes that have at least one overlapping gene
load("data/general/overlappingGenes.rda")

neuroNonOveGen <- geneNames
## remove the non automosal and the patch genes, don't have the coverage on them
neuroNonOveGen <- neuroNonOveGen[-c(147:158),]

# neuroNonOveGen <- geneNames[-which(geneNames[,1] %in% rownames(as.data.frame(listNonOve))),]
# neuroNonOveGen <- neuroNonOveGen[order(neuroNonOveGen$chromosome_name),]

novelTransRegion(neuroNonOveGen[3,],ensembl,10,"PUTM")



library(doParallel)
library(foreach)
library(GenomicRanges)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(5)
clusterExport(cl, c("novelTransRegion","getBM","subsetByOverlaps"))
registerDoParallel(cl)
getDoParWorkers()

## remove gene that comes from the patch HG987_PATCH - gene KCNJ12
## neuroNonOveGen <- neuroNonOveGen[c(-58),]


# neuroNonOveGen <- neuroNonOveGen[-c(58,59,60,61,62),]



start <- Sys.time()
novelRegions <- foreach(i=1:nrow(neuroNonOveGen),.combine=rbind,.verbose=F)%dopar%novelTransRegion(neuroNonOveGen[i,],ensembl,10,"PUTM")
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)


neuroGenes <- getBM(attributes=c("ensembl_gene_id","external_gene_id","gene_biotype","source","status"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=unlist(lapply(novelRegions,function(x){return(x$gene)})), mart=ensembl)


rownames(novelRegions) <- novelRegions[,1]
novelRegions <- novelRegions[,-1]
rownames(neuroGenes) <- neuroGenes[,1]
neuroGenes <- neuroGenes[,-1]
neuroGenes <- neuroGenes[rownames(novelRegions),]
identical(rownames(neuroNonOveGen),rownames(neuroGenes))
neuroGenes <- cbind(novelRegions,neuroGenes)

neuroGenes$overlapGene <- FALSE
neuroGenes$overlapGene[which(rownames(neuroGenes) %in% rownames(as.data.frame(listNonOve)))] <- TRUE
save(neuroGenes,file="data/results/novelIntragenicRegions.PUTM.rda")

### SNIG
rm(neuroGenes,novelRegions)

start <- Sys.time()
novelRegions <- foreach(i=1:nrow(neuroNonOveGen),.combine=rbind,.verbose=F)%dopar%novelTransRegion(neuroNonOveGen[i,],ensembl,10,tissue="SNIG")
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)

neuroGenes <- getBM(attributes=c("ensembl_gene_id","external_gene_id","gene_biotype","source","status"),
                    verbose = T,
                    filters="ensembl_gene_id",
                    values=novelRegions[,1], mart=ensembl)

rownames(novelRegions) <- novelRegions[,1]
novelRegions <- novelRegions[,-1]
rownames(neuroGenes) <- neuroGenes[,1]
neuroGenes <- neuroGenes[,-1]
neuroGenes <- neuroGenes[rownames(novelRegions),]
identical(rownames(neuroNonOveGen),rownames(neuroGenes))
neuroGenes <- cbind(novelRegions,neuroGenes)

neuroGenes$overlapGene <- FALSE
neuroGenes$overlapGene[which(rownames(neuroGenes) %in% rownames(as.data.frame(listNonOve)))] <- TRUE
save(neuroGenes,file="data/results/novelIntragenicRegions.SNIG.rda")



neuroGenes[which(neuroGenes$external_gene_id %in% "APP"),]
## we now try to see whether we could predict new regions checking in a recent version of ensembl


cl <- makeCluster(5)
clusterExport(cl, c("novelTransRegion","getBM","subsetByOverlaps"))
registerDoParallel(cl)
getDoParWorkers()

ensemblv75 <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Feb2014.archive.ensembl.org",
                      dataset="hsapiens_gene_ensembl")

start <- Sys.time()
novelRegionsv75 <- foreach(i=1:nrow(neuroNonOveGen),.combine=rbind,.verbose=F)%dopar%novelTransRegion(neuroNonOveGen[i,],ensemblv75)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start

sum(unlist(novelRegionsv75[,3]),na.rm=T)
stopCluster(cl)
rm(cl)




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

load_all()

load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
IDs=PUTM$A.CEL_file


plotReadDepth(gene="ENSG00000186868",ensembl=ensembl,IDs=IDs)

load("data/general/sampleInfo.rda")
SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified =="SNIG"),]
IDs=SNIG$A.CEL_file

plotReadDepth(gene="ENSG00000186868",ensembl=ensembl,IDs=IDs)

head(novelRegions,20)






### we now try to identify the actual regions:
load(file="data/results/novelIntragenicRegions.SNIG.rda")

head(neuroGenes)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")
tissue <- "SNIG"


##select only the genes that have potentially novel intragenic transcribed regions
neuroGenesTmp <- neuroGenes[which(neuroGenes$coveredNoAnn.NA>0),]

library(doParallel)
library(foreach)
library(GenomicRanges)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(1)
clusterExport(cl, c("getUnAnnotatedRegions","getBM","subsetByOverlaps","GRanges","countOverlaps"))
registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
novelIntraRegion <- foreach(i=1:nrow(neuroGenesTmp),.verbose=F)%dopar%getUnAnnotatedRegions(rownames(neuroGenesTmp)[i],ensembl,"SNIG")
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)


novelIntraRegion <- GRangesList(novelIntraRegion)
names(novelIntraRegion) <- rownames(novelIntraRegion)

save(novelIntraRegion,file="data/results/novelIntragenicRegionsLoc.SNIG.rda")
getwd()

