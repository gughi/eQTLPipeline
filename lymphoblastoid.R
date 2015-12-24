## Comaparison Lymphoblastoid

## We removed genes and SNP that are not present in our eQTL analysi

cmd <- "cut -f1,3 /home/seb/projectsR/eQTLPipeline/data/general/Lappailanen/EUR373.gene.cis.FDR5.all.rs137.txt | cut -d '.' -f1 | cut -f1,2"
lym.eQTLs.pairs <- read.delim(pipe(cmd),header=T)
lym.eQTLs <- read.delim(file="data/general/Lappailanen/EUR373.gene.cis.FDR5.all.rs137.txt",header=T)
## we get the slope | rvalue  : Spearman rank correlation rho (calculated from linear regression slope). The sign denotes the direction of the nonreference allele (i.e. rvalue<0 means that nonreference allele has lower expression) 
lym.eQTLs.pairs$slope <- lym.eQTLs$rvalue 
## nrow(lym.eQTLs)
## [1] 419983

length(unique(lym.eQTLs.pairs$GENE_ID))
## 3259

load("data/general/Lappailanen/filteredFormattedRSnumbers.rda")

rsPos <- rsPos[-which(duplicated(rsPos$rsNumber)),]
rownames(rsPos) <- rsPos$rsNumber 
rsPos$rsNumber <- NULL
head(rsPos)
head(rsPos[as.character(lym.eQTLs.pairs$SNP_ID),])
lym.eQTLs.pairs$chrPos <- rsPos[as.character(lym.eQTLs.pairs$SNP_ID),"formattedChrPos"]


rm(cmd)
cmd <- "cut -f1 -d ' ' /home/adai/genotyped/imputed_v3/imputed.info"
snpsBrain <- read.delim(pipe(cmd))
rm(cmd)

lym.eQTLs.pairs <- lym.eQTLs.pairs[which(lym.eQTLs.pairs$chrPos %in% snpsBrain[,1]),]


load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")


lym.eQTLs.pairs <- lym.eQTLs.pairs[which(lym.eQTLs.pairs$GENE_ID %in% as.character(colnames(resids))),]

load("data/results/finaleQTLs/geneExonic.unsentinalised.rda")

eQTLPUTMExun <- eQTLPUTMExun[which(eQTLPUTMExun[,7] < 0.05),]      
eQTLSNIGExun <- eQTLSNIGExun[which(eQTLSNIGExun[,7] < 0.05),]     

## remove genes that are not present in both studies
lymGene <- read.delim(pipe("cut -f1 data/general/Lappailanen/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt | cut -d '.' -f1"))

eQTLPUTMExun <- eQTLPUTMExun[which(eQTLPUTMExun[,2] %in% lymGene[,1]),]      

chrPos.chr <- gsub("chr","",sapply(strsplit(as.character(eQTLPUTMExun[,1]), split=":",fixed=T ), function(x) x[1] ))
chrPos.pos <- sapply(strsplit(as.character(eQTLPUTMExun[,1]), split=":",fixed=T ), function(x) x[2] )  

eQTLGeneExonic <- data.frame(paste(chrPos.chr,chrPos.pos,sep=":"))
eQTLGeneExonic$genes <- eQTLPUTMExun[,2]
eQTLGeneExonic$beta <- eQTLPUTMExun[,6]  
colnames(eQTLGeneExonic)[1] <- "snps"

rm(chrPos.chr,chrPos.pos,cmd)
rsPos[as.character(head(lym.eQTLs.pairs$SNP_ID)),]

lym.eQTLs.pairs$chrPos <- rsPos[as.character(lym.eQTLs.pairs$SNP_ID),"formattedChrPos"]
head(lym.eQTLs.pairs)

length(intersect(paste0(eQTLGeneExonic$snps,eQTLGeneExonic$genes),paste0(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID)))

tmp <- intersect(paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),paste(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID,sep="_"))

length(unique(sapply(strsplit(as.character(tmp), split="_",fixed=T ), function(x) x[2] )))
rm(tmp)

library(gplots)
venn(list(Putamen=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),Lymphoblastoid=paste(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID,sep="_")))

venn(list(Putamen=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sign(as.numeric(eQTLGeneExonic$beta)),sep="_"),
          Lymphoblastoid=paste(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID,sign(as.numeric(lym.eQTLs.pairs$slope)),sep="_")))

## check the intersectio of genes 

tmp <- intersect(paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sign(as.numeric(eQTLGeneExonic$beta)),sep="_"),
        paste(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID,sign(as.numeric(lym.eQTLs.pairs$slope)),sep="_"))

length(unique(sapply(strsplit(as.character(tmp), split="_",fixed=T ), function(x) x[2] )))




eQTLSNIGExun <- eQTLSNIGExun[which(eQTLSNIGExun[,2] %in% lymGene[,1]),]      

chrPos.chr <- gsub("chr","",sapply(strsplit(as.character(eQTLSNIGExun[,1]), split=":",fixed=T ), function(x) x[1] ))
chrPos.pos <- sapply(strsplit(as.character(eQTLSNIGExun[,1]), split=":",fixed=T ), function(x) x[2] )  

eQTLGeneExonic <- data.frame(paste(chrPos.chr,chrPos.pos,sep=":"))
eQTLGeneExonic$genes <- eQTLSNIGExun[,2]
eQTLGeneExonic$beta <- eQTLSNIGExun[,6]  
colnames(eQTLGeneExonic)[1] <- "snps"

rm(chrPos.chr,chrPos.pos,cmd)
rsPos[as.character(head(lym.eQTLs.pairs$SNP_ID)),]

lym.eQTLs.pairs$chrPos <- rsPos[as.character(lym.eQTLs.pairs$SNP_ID),"formattedChrPos"]
head(lym.eQTLs.pairs)

length(intersect(paste0(eQTLGeneExonic$snps,eQTLGeneExonic$genes),paste0(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID)))

tmp <- intersect(paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),paste(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID,sep="_"))

length(unique(sapply(strsplit(as.character(tmp), split="_",fixed=T ), function(x) x[2] )))


venn(list(Putamen=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),Lymphoblastoid=paste(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID,sep="_")))

venn(list(SubstantiaNigra=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sign(as.numeric(eQTLGeneExonic$beta)),sep="_"),
          Lymphoblastoid=paste(lym.eQTLs.pairs$chrPos,lym.eQTLs.pairs$GENE_ID,sign(as.numeric(lym.eQTLs.pairs$slope)),sep="_")))






load("data/results/finaleQTLs/geneExonic.unsentinalised.rda")

eQTLPUTMExun <- eQTLPUTMExun[which(eQTLPUTMExun[,7] < 0.05),]      
eQTLSNIGExun <- eQTLSNIGExun[which(eQTLSNIGExun[,7] < 0.05),]     

chrPos.chr <- gsub("chr","",sapply(strsplit(as.character(eQTLPUTMExun[,1]), split=":",fixed=T ), function(x) x[1] ))
chrPos.pos <- sapply(strsplit(as.character(eQTLPUTMExun[,1]), split=":",fixed=T ), function(x) x[2] )  
eQTLGeneExonic <- data.frame(paste(chrPos.chr,chrPos.pos,sep=":"))
eQTLGeneExonic$genes <- eQTLPUTMExun[,2]
eQTLGeneExonic$beta <- eQTLPUTMExun[,6]  
colnames(eQTLGeneExonic)[1] <- "snps"
rm(chrPos.chr,chrPos.pos)

GTEx <- read.delim("data/general/GTEx/Brain_Putamen_basal_ganglia_Analysis.snpgenes")

gene <- unlist(sapply(strsplit(as.character(GTEx$gene), split=".",fixed=T ), function(x) x[1] ))

GTEx$gene <- gene
GTEx$snp <- paste(GTEx$snp_chrom,GTEx$snp_pos,sep=":")



load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")

GTEx <- GTEx[which(GTEx$gene %in% as.character(colnames(resids))),]


cmd <- "cut -f1 -d ' ' /home/adai/genotyped/imputed_v3/imputed.info"
snpsBrain <- read.delim(pipe(cmd))
rm(cmd)

GTEx <- GTEx[which(GTEx$snp %in% snpsBrain[,1]),]


venn(list(UKBEC=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),
          GTEx=paste(GTEx$snp,GTEx$gene,sep="_")))


cmd <- "cut -f1 /home/seb/projectsR/eQTLPipeline/data/general/GTEx/Brain_Putamen_basal_ganglia_Analysis.expr.txt | cut -d '.' -f1"
GTExGenes <- read.delim(pipe(cmd),header=T)


eQTLGeneExonic <- eQTLGeneExonic[which(eQTLGeneExonic$genes %in% GTExGenes$Id),]

venn(list(UKBEC=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),
          GTEx=paste(GTEx$snp,GTEx$gene,sep="_")))

venn(list(UKBEC=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sign(as.numeric(eQTLGeneExonic$beta)),sep="_"),
          GTEx=paste(GTEx$snp,GTEx$gene,sign(as.numeric(GTEx$beta)),sep="_")))

tmp <- intersect(paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),paste(GTEx$snp,GTEx$gene,sep="_"))

length(unique(sapply(strsplit(as.character(tmp), split="_",fixed=T ), function(x) x[2] )))







load("data/results/finaleQTLs/geneExonic.unsentinalised.rda")

eQTLPUTMExun <- eQTLPUTMExun[which(eQTLPUTMExun[,7] < 0.05),]      
eQTLSNIGExun <- eQTLSNIGExun[which(eQTLSNIGExun[,7] < 0.05),]     

chrPos.chr <- gsub("chr","",sapply(strsplit(as.character(eQTLPUTMExun[,1]), split=":",fixed=T ), function(x) x[1] ))
chrPos.pos <- sapply(strsplit(as.character(eQTLPUTMExun[,1]), split=":",fixed=T ), function(x) x[2] )  
eQTLGeneExonic <- data.frame(paste(chrPos.chr,chrPos.pos,sep=":"))
eQTLGeneExonic$genes <- eQTLPUTMExun[,2]
eQTLGeneExonic$beta <- eQTLPUTMExun[,6]  
colnames(eQTLGeneExonic)[1] <- "snps"
rm(chrPos.chr,chrPos.pos)

GTEx <- read.delim("data/general/GTEx/Brain_Hippocampus_Analysis.snpgenes")

gene <- unlist(sapply(strsplit(as.character(GTEx$gene), split=".",fixed=T ), function(x) x[1] ))

GTEx$gene <- gene
GTEx$snp <- paste(GTEx$snp_chrom,GTEx$snp_pos,sep=":")

load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")

GTEx <- GTEx[which(GTEx$gene %in% as.character(colnames(resids))),]


# cmd <- "cut -f1 -d ' ' /home/adai/genotyped/imputed_v3/imputed.info"
cmd <- "cut -f1 -d ' ' /Volumes/NONAME/imputed/imputed.dosage"
snpsBrain <- read.delim(pipe(cmd))
rm(cmd)

GTEx <- GTEx[which(GTEx$snp %in% snpsBrain[,1]),]

library(gplots)
venn(list(UKBEC=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),
          GTEx=paste(GTEx$snp,GTEx$gene,sep="_")))


cmd <- "cut -f1 /home/seb/projectsR/eQTLPipeline/data/general/GTEx/Brain_Putamen_basal_ganglia_Analysis.expr.txt | cut -d '.' -f1"
GTExGenes <- read.delim(pipe(cmd),header=T)


eQTLGeneExonic <- eQTLGeneExonic[which(eQTLGeneExonic$genes %in% GTExGenes$Id),]

venn(list(UKBEC=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),
          GTEx=paste(GTEx$snp,GTEx$gene,sep="_")))

venn(list(UKBEC=paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sign(as.numeric(eQTLGeneExonic$beta)),sep="_"),
          GTEx=paste(GTEx$snp,GTEx$gene,sign(as.numeric(GTEx$beta)),sep="_")))

tmp <- intersect(paste(eQTLGeneExonic$snps,eQTLGeneExonic$genes,sep="_"),paste(GTEx$snp,GTEx$gene,sep="_"))

length(unique(sapply(strsplit(as.character(tmp), split="_",fixed=T ), function(x) x[2] )))















