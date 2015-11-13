## For Jana
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


tmp <- unlist(lapply(strsplit(as.character(mapExon$geneID),"_")
                            ,function(x){c(x[1])}))

mapExonNeuro <- mapExon[which(tmp %in% geneNames[,1]),]

# potExExJun <- expand.grid(mapExonNeuro$exonID,mapExonNeuro$exonID)
# potExExJun <- paste(potExExJun[,1],potExExJun[,2],sep = "_")

idx <- unique(which(expr$Exon1ID %in% mapExonNeuro$exonID),which(expr$Exon2ID %in% mapExonNeuro$exonID))
exprExExJun <- expr[idx,]

mapExExNeuroGenes <- mapExon[which(as.character(mapExon$exonID) %in% unique(c(as.character(exprExExJun$Exon1ID),as.character(exprExExJun$Exon2ID)))),]
exprExExNeuroGenes <- exprExExJun
save(exprExExNeuroGenes,mapExExNeuroGenes,file="data/expr/rawCounts/genic/exExJunNeuroGenes.rda")


intersect(as.character(mapExon$exonID),as.character(exprExExNeuroGenes$Exon1ID))
trans1 <- mapExon[which(as.character(mapExon$exonID) %in% as.character(exprExExNeuroGenes$Exon1ID)),c("exonID","transID")]
trans2 <- mapExExNeuroGenes[which(as.character(mapExExNeuroGenes$exonID) %in% as.character(exprExExNeuroGenes$Exon2ID)),c("exonID","transID")]


trans1 <- trans1[as.character(exprExExNeuroGenes$Exon1ID),"transID"]
trans2 <- trans2[as.character(exprExExNeuroGenes$Exon2ID),"transID"]

trans1 <- unlist(lapply(strsplit(as.character(trans1),"_")
                     ,function(x){c(x[2])}))

trans2 <- unlist(lapply(strsplit(as.character(trans2),"_")
                        ,function(x){c(x[2])}))

exprExExNeuroGenes$potNove <- trans1==trans2
save(exprExExNeuroGenes,mapExExNeuroGenes,file="data/expr/rawCounts/genic/exExJunNeuroGenes.rda")



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


