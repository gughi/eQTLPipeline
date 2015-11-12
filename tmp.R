## For Jana
load("data/expr/rawCounts/genic/fullExExJun.rda")
gene <- read.delim("data/general/NP4.6c.raw.n.bed")
head(mapExon)





## I obtain the gene list
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



head(geneNames)
## remove all the LRG genes
geneNames <- geneNames[which(geneNames$gene_biotype %in% "protein_coding"),]  

head(geneNames)

library(GenomicRanges)

GR <- GRanges(seqnames = Rle(geneNames$chromosome_name),
              ranges = IRanges(start=geneNames$start_position,
                               end = geneNames$end_position,
                               names = geneNames$Ensembl.Gene.ID))










