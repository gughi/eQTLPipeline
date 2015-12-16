## For Rohan Locus
 
## Rohan is interested in the locus is MOBP and the hg19 co-ordinates: chr3:39,509,064-39,570,988.
## eQTL reported in CRBL for 3:39523003 

## The gene in ensembl v72 is
## 3:39508689-39570970:1
MOBP <-  "ENSG00000168314"

## this is what we are seeing for the sentinal SNPs
eQTLPUTM <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")

genes <- unlist(lapply(strsplit(as.character(eQTLPUTM$gene),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))


eQTLPUTM[which(genes %in% MOBP),]

## SNIG
eQTLSNIG <- read.delim("data/results/finaleQTLs/exons.SNIG.txt",sep=" ")

genes <- unlist(lapply(strsplit(as.character(eQTLSNIG$gene),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))


eQTLSNIG[which(genes %in% MOBP),]

# snps                 gene    tstat       pvalue        FDR      beta      myFDR degree
# 710  chr3:39510573 ENSG00000168314:E030 4.794916 1.113397e-05 0.05082659 0.3356933 0.05082659      1
# 2104 chr3:39510573 ENSG00000168314:E019 4.801486 1.087344e-05 0.04963723 0.4622903 0.04963723      1
# 2707 chr3:39510573 ENSG00000168314:E031 4.995177 5.382107e-06 0.02456932 0.3292004 0.02456932      1
# 4298 chr3:39510573 ENSG00000168314:E032 5.144414 3.110221e-06 0.01419816 0.3274482 0.01419816      1
# 4301 chr3:39510573 ENSG00000168314:E017 4.633457 1.984387e-05 0.09058726 0.4109740 0.09058726      1
# 4540 chr3:39510573 ENSG00000168314:E018 5.114629 3.471456e-06 0.01584720 0.4714849 0.01584720      1

## rs28764178, rs number for the SNP we found eQTL
 

## So for the SNP that rohan is interested and the SNP we are finding eQTLs there is a distance only of 12430
## 39523003-39510573

## annotation for the exonss:
transcriptomeInfo <- read.csv("data/general/transcriptomeInfo.csv",row.names=4)
head(transcriptomeInfo)
transcriptomeInfo[as.character(eQTLSNIG[which(genes %in% MOBP),2]),]

# ENSG00000168314:E030 39565935 39565985    51
# ENSG00000168314:E019 39544372 39544405    34
# ENSG00000168314:E031 39565986 39566047    62
# ENSG00000168314:E032 39566048 39566373   326
# ENSG00000168314:E017 39544026 39544367   342
# ENSG00000168314:E018 39544368 39544371     4



rm(eQTLPUTM,eQTLSNIG,genes,transcriptomeInfo)

load("data/results/finaleQTLs/exons.unsentinalised.formatted.SNIG.txt")

genes <- unlist(lapply(strsplit(as.character(eQTLSNIG[,2]),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

unique(eQTLSNIG[which(genes %in% MOBP),1])

rm(eQTLSNIG,genes)

## Exon exon junctions

load("data/results/finaleQTLs/eQTL.ExExJun.PUTM.rda")
eQTL.PUTM[which(eQTL.PUTM$geneID %in% MOBP),]

load("data/results/finaleQTLs/eQTL.ExExJun.SNIG.rda")
eQTL.SNIG[which(eQTL.SNIG$geneID %in% MOBP),]

# snps        exExID    tstat       pvalue      FDR      beta    myFDR degree          geneID                         transID chrExon1
# 373351_373355 chr3:39510573 373351_373355 5.034244 4.664842e-06 0.021295 0.4006396 0.021295      1 ENSG00000168314 ENST00000424090;ENST00000479860        3
# startExon1 endExon1 chrExon2 startExon2 endExon2 differentGene distanceExons numOfTrans hasUniExo external_gene_id   gene_biotype
# 373351_373355   39565733 39565817        3   39565934 39566047         FALSE           117          1     FALSE             MOBP protein_coding

plotLoceQTLs(gene=MOBP,gen="hg19",ensembl=ensembl,IDs=IDs)


rm(eQTL.PUTM,eQTL.SNIG)

load_all()

load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
IDs <- gsub("/","_",PUTM$U.SD_No)
IDs=PUTM$A.CEL_file

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

plotReadDepth(gene=MOBP,gen="hg19",ensembl=ensembl,IDs=IDs)




load(paste0("/home/seb/eQTL/snps/byGene/",MOBP,".rda"))
snp <- "chr3:39510573"
markers <- markers[snp,]

load("data/general/sampleInfo.rda")
SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified =="SNIG"),]
IDs <- gsub("/","_",SNIG$U.SD_No)
tmp <- markers[,as.character(IDs)]
names(tmp) <- as.character(SNIG$A.CEL_file)
markers <- list(info=markers[,c(1:6)],genotype=tmp)
rm(IDs,tmp)

table(round(as.numeric(markers$genotype)))
genotype=markers
IDs=SNIG$A.CEL_file


library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

plotReadDepth(gene=MOBP,gen="hg19",ensembl=ensembl,IDs=IDs)

load("data/results/finaleQTLs/eQTL.ExExJun.SNIG.rda")
eQTL.SNIG[which(eQTL.SNIG$geneID %in% MOBP),]

plotLoceQTLs(gene = MOBP,ensembl = ensembl,IDs = IDs,genotype = genotype,
             highLight=eQTL.SNIG[intersect(which(eQTL.SNIG$snps %in% snp),which(eQTL.SNIG$geneID %in% MOBP)),
                                 c("chrExon1","startExon1","endExon1","startExon2","endExon2","exExID")])
  















