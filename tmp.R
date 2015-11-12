## For Jana
load("data/expr/rawCounts/genic/fullExExJun.rda")
gene <- read.delim("data/general/NP4.6c.raw.n.bed")
head(mapExon)





## I obtain the gene list
gene <- unique(gene$LOCUS)


library(biomaRt)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


<<<<<<< HEAD

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








=======
# ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Feb2014.archive.ensembl.org",
#                       dataset="hsapiens_gene_ensembl")

plotLoceQTLs(gene = gene,ensembl = ensembl,IDs = IDs,genotype = genotype)


## an example of exon-exon eQTL
eQTL.PUTM[1,]

gene <- "ENSG00000230051"
snp <- "chr22:27731294"
##

## load(paste0("data/snps/byGene/",gene,".rda"))
##
load(paste0("/home/seb/eQTL/snps/byGene/",gene,".rda"))
markers <- markers[snp,]

# cooment 


load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
IDs <- gsub("/","_",PUTM$U.SD_No)
tmp <- markers[,as.character(IDs)]
names(tmp) <- as.character(PUTM$A.CEL_file)
markers <- list(info=markers[,c(1:6)],genotype=tmp)
rm(IDs,tmp)

table(round(as.numeric(markers$genotype)))
genotype=markers
IDs=PUTM$A.CEL_file

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


# ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Feb2014.archive.ensembl.org",
#                       dataset="hsapiens_gene_ensembl")

plotLoceQTLs(gene = gene,ensembl = ensembl,IDs = IDs,genotype = genotype)



#### 
gene <- "ENSG00000163930"
snp <- "chr3:53101224"
##

## load(paste0("data/snps/byGene/",gene,".rda"))
##
load(paste0("/home/seb/eQTL/snps/byGene/",gene,".rda"))
markers <- markers[snp,]

# cooment 


load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
IDs <- gsub("/","_",PUTM$U.SD_No)
tmp <- markers[,as.character(IDs)]
names(tmp) <- as.character(PUTM$A.CEL_file)
markers <- list(info=markers[,c(1:6)],genotype=tmp)
rm(IDs,tmp)

table(round(as.numeric(markers$genotype)))
genotype=markers
IDs=PUTM$A.CEL_file

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


# ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Feb2014.archive.ensembl.org",
#                       dataset="hsapiens_gene_ensembl")

plotLoceQTLs(gene = gene,ensembl = ensembl,IDs = IDs,genotype = genotype)


>>>>>>> 86eddcddb4a480eaebdf1b173ca9fbb2b31b0965


