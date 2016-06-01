#' ---
#' title: "Enrichment overlap cis-eQTLs intergenic transcribed regions and GWAS hits"
#' author: "Manuel Sebastian Guelfi"
#' ---
#'
#' The report shows the overlap between GWAS and the cis-eQTL intergenic region and exons, we compare with exons
#' because we think that this is the more similar expression
#' 
#'      

#+ echo=FALSE
suppressWarnings(library(knitr))
opts_chunk$set(echo=FALSE)
## need also to set the directory for the entire script



#' Downloaded only three histone marks as suggested by Sarah G. She said these three are the most representative histone modifciation. 
#' 1. H3K4Me1- mark of enhancer (active enhancer when with H3K27Ac)
#' 2. H3K4Me3- mark of promoter
#' 3. H3K27Ac

#' Downloaded from:
#' http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/gappedPeak/
#' on the 9th of May 2016
suppressMessages(library(devtools))
suppressMessages(library(GenomicRanges))
setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()
library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
opts_knit$set(root.dir=dirname(path$networkPathname))
rm(path)


H3K4me1 <- read.delim("data/Roadmap/E074-H3K4me1.imputed.gappedPeak.bed.gPk",header=F)
H3K4me3 <- read.delim("data/Roadmap/E074-H3K4me3.imputed.gappedPeak.bed.gPk",header=F)
H3K27ac <- read.delim("data/Roadmap/E074-H3K27ac.imputed.gappedPeak.bed.gPk",header=F)
DNase <- read.delim("data/Roadmap/E074-DNase.imputed.gappedPeak.bed.gPk",header=F)

#' convert all the ranges in Genomics range object

DNase <- GRanges(DNase[,1], IRanges(DNase[,2], DNase[,3]))
sum(width(DNase))
#'  Proportion of genoic regions for DNase: 231,647,494
H3K4me1 <- GRanges(H3K4me1[,1], IRanges(H3K4me1[,2], H3K4me1[,3]))
sum(width(H3K4me1))
#'  Proportion of genoic regions for H3K4me1: 436,818,199
H3K4me3 <- GRanges(H3K4me3[,1], IRanges(H3K4me3[,2], H3K4me3[,3]))
sum(width(H3K4me3))
#'  Proportion of genoic regions for H3K4me3: 101,058,388
H3K27ac <- GRanges(H3K27ac[,1], IRanges(H3K27ac[,2], H3K27ac[,3]))
sum(width(H3K27ac))
#'  Proportion of genoic regions for H3K27ac: 238,051,848

nigra <- list(DNase=DNase,H3K27ac=H3K27ac,H3K4me3=H3K4me3,H3K4me1=H3K4me1)
rm(H3K27ac,H3K4me3,H3K4me1,DNase)

#' So first thing I checked how much was the overlap, for example between two diffent histone marks for methylation. 
over <- overlapsAny(nigra$H3K4me3,nigra$H3K4me1,ignore.strand=T )
table(over)
#' FALSE  TRUE 
#' 1641 37597 
#' Most of the regions overlap and I was concern. Discuss it with Jana and the reasons can be multiple:
#' - We are measuring different types of cells, it's human brain tissue. So some of the histones 
#' overlap because they are measuring the "peaks" for different cells
#' - Since this histones can physically very close and that will get silmilar peaks when two peaks are 
#' in contact with each other
#' - The histone modifications peaks are the product of different experiments
#' In conclusion, we should not pay to much attention at the fact that they are overlapping each other.
rm(over)

load("data/results/novelIntergenicRegionsFinal/intergenicRegion.Ann.rda")

eQTLsRanges <- GRanges(eQTLs.SNIG$chr, IRanges(eQTLs.SNIG$start, eQTLs.SNIG$end))

# we get overlap for each region
over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})


tmp <- cbind(ann=eQTLs.SNIG$ourAnn,
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})





# we get the summary of the overlap
summOver <- lapply(over, function(x){
                          tmp <- cbind(ann=eQTLs.SNIG$ourAnn,marker=x)
                           aggregate(marker ~ ann,data=tmp,table) 
                         })
rm(tmp)


sapply(seq_along(summOver), function(i){
  barplot(t(summOver[[i]][,2]),beside = T,col = 1:2,main=names(summOver[i]),
          legend.text = c("Non-overlapping","Overlapping"),args.legend = list(x = "topleft"),
          names.arg = summOver[[i]][,1])})


sapply(seq_along(summOver), function(i){
  barplot(summOver[[i]][,2][,2]/apply(summOver[[i]][,2],1,sum),beside = T,col = 1:4,main=names(summOver[i]),
          names.arg = summOver[[i]][,1],ylab = "ratio")})

tmp1 <- as.matrix(as.data.frame(lapply(summOver,function(x){x[,2][,"TRUE"]/apply(x[,2],1,sum)})))
rownames(tmp1) <- summOver[[1]][,1]
rm(summOver)



#' Next check is get how much is the overlap between the expressed annotated coding regions (exons) and check how much is the overlap

load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.SNIG")
rm(covs,PUTM,SNIG)
genes <- rownames(RPKM.cqn)
rm(RPKM.cqn)
exonDef <- read.csv("data/general/transcriptomeInfo.csv")
exonDef <- exonDef[which(as.character(genes) %in% as.character(exonDef$names)),]
rm(genes)

genes <- unlist(lapply(strsplit(as.character(exonDef$names),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))
geneDef <- read.delim("data/general/ensemblGenes.txt")

exonDef <- cbind(exonDef,Chromosome.Name=geneDef[match(as.character(genes),as.character(geneDef$Ensembl.Gene.ID)),"Chromosome.Name"])

library("biomaRt")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","gene_biotype"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=as.character(genes), mart=ensembl)

#head(geneNames)
exonDef <- cbind(exonDef,biotype=geneNames[match(as.character(genes),as.character(geneNames$ensembl_gene_id)),"gene_biotype"])
bp <- barplot(table(exonDef$biotype),main="Gene biotype for the brain expressed exons - SNIG")
text(bp , 0,labels = table(exonDef$biotype),cex=1,pos=3)

rm(genes,bp,geneDef)

## we create the ranges for the exons
exonRanges <- GRanges(paste0("chr",exonDef$Chromosome.Name), IRanges(exonDef$start, exonDef$end))

rm(over)
over <- lapply(nigra, function(x){overlapsAny(exonRanges,x)})

tmp <- cbind(biotype=as.character(exonDef$biotype),
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})


## get the overlap based on the biotype

tmp <- cbind(tmp,allMarks=over$allMarks)

summOver <- list(DNAse=aggregate(DNase ~ biotype,data=tmp,table),
       H3K27ac=aggregate(H3K27ac ~ biotype,data=tmp,table),
       H3K4me1=aggregate(H3K4me1 ~ biotype,data=tmp,table),
       H3K4me3=aggregate(H3K4me3 ~ biotype,data=tmp,table),
       allMarks=aggregate(allMarks ~ biotype,data=tmp,table))

#' We now visualise the overlap for the brain annotated expressed exons divided by biotype 
barplot(sapply(seq_along(summOver), function(i){summOver[[i]][,2][,"TRUE"]/
    apply(summOver[[i]][,2],1,sum)}),
        names.arg = names(summOver),ylab = "ratio",
    main="Overlap annotated expressed exons by biotype - SNIG",beside = T,
    legend.text = summOver[[1]][,1],col=1:6,ylim=c(0,0.7),args.legend = list(x = "topleft"))


#' We are not convince by the overlap of the lincRNA, specially with the histone mark H3K4me3 that should be 
#' a marker for lincRNA. So check whether increasing the range for the overlap can improve the overlap

over <- lapply(nigra, function(x){overlapsAny(exonRanges,x,maxgap = 5000)})

tmp <- cbind(biotype=as.character(exonDef$biotype),
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)

over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})


## get the overlap based on the biotype

tmp <- cbind(tmp,allMarks=over$allMarks)

summOver <- list(DNAse=aggregate(DNase ~ biotype,data=tmp,table),
                 H3K27ac=aggregate(H3K27ac ~ biotype,data=tmp,table),
                 H3K4me1=aggregate(H3K4me1 ~ biotype,data=tmp,table),
                 H3K4me3=aggregate(H3K4me3 ~ biotype,data=tmp,table),
                 allMarks=aggregate(allMarks ~ biotype,data=tmp,table))

#' We now visualise the overlap for the brain annotated expressed exons divided by biotype 
barplot(sapply(seq_along(summOver), function(i){summOver[[i]][,2][,"TRUE"]/
    apply(summOver[[i]][,2],1,sum)}),
    names.arg = names(summOver),ylab = "ratio",
    main="Overlap annotated expressed exons by biotype - SNIG",beside = T,
    legend.text = summOver[[1]][,1],col=1:6,ylim=c(0,0.7),args.legend = list(x = "topleft"))


#' using a 5Kb as a distance gap would not made any difference in the overlap.


rm(summOver)

summOver <- lapply(over,table)

barplot(sapply(seq_along(summOver), function(i){summOver[[i]]["TRUE"]/sum(summOver[[i]])}),
        names.arg = names(summOver),ylab = "ratio", main="Overlap with expressed exons SNIG")



tmp1 <- as.matrix(rbind(tmp1,annExons=as.data.frame(lapply(summOver, function(x){x["TRUE"]/sum(x)}))[1,]))

barplot(tmp1,beside = T, 
        args.legend = list(x = "topleft"),legend.text = rownames(tmp1),col=1:5,ylab = "ratio",main="Overlap SNIG")





#' It seems we have low overlap specially for independent regios and novel exons.


#'We apply the same thing to Putamen. We use **Caudate** as most similar to **Putamen**

rm(list=ls())

H3K4me1 <- read.delim("data/Roadmap/E068-H3K4me1.imputed.gappedPeak.bed.gPk",header=F)
H3K4me3 <- read.delim("data/Roadmap/E068-H3K4me3.imputed.gappedPeak.bed.gPk",header=F)
H3K27ac <- read.delim("data/Roadmap/E068-H3K27ac.imputed.gappedPeak.bed.gPk",header=F)
DNase <- read.delim("data/Roadmap/E068-DNase.imputed.gappedPeak.bed.gPk",header=F)

#' convert all the ranges in Genomics range object

DNase <- GRanges(DNase[,1], IRanges(DNase[,2], DNase[,3]))
sum(width(DNase))
#'  Proportion of genoic regions for DNase: 245,882,474
H3K4me1 <- GRanges(H3K4me1[,1], IRanges(H3K4me1[,2], H3K4me1[,3]))
sum(width(H3K4me1))
#'  Proportion of genoic regions for H3K4me1: 434,696,531
H3K4me3 <- GRanges(H3K4me3[,1], IRanges(H3K4me3[,2], H3K4me3[,3]))
sum(width(H3K4me3))
#'  Proportion of genoic regions for H3K4me3: 98,297,179
H3K27ac <- GRanges(H3K27ac[,1], IRanges(H3K27ac[,2], H3K27ac[,3]))
sum(width(H3K27ac))
#'  Proportion of genoic regions for H3K27ac: 245,093,759

caudate <- list(DNase=DNase,H3K27ac=H3K27ac,H3K4me3=H3K4me3,H3K4me1=H3K4me1)
rm(H3K27ac,H3K4me3,H3K4me1,DNase)

#' So first thing I checked how much was the overlap, for example between two diffent histone marks for methylation. 
over <- overlapsAny(caudate$H3K4me3,caudate$H3K4me1,ignore.strand=T )
table(over)
#' FALSE  TRUE 
#' 1552 36027 
#' Most of the regions overlap and I was concern. Discuss it with Jana and the reasons can be multiple:
#' - We are measuring different types of cells, it's human brain tissue. So some of the histones 
#' overlap because they are measuring the "peaks" for different cells
#' - Since this histones can physically very close and that will get silmilar peaks when two peaks are 
#' in contact with each other
#' - The histone modifications peaks are the product of different experiments
#' In conclusion, we should not pay to much attention at the fact that they are overlapping each other.
rm(over)

load("data/results/novelIntergenicRegionsFinal/intergenicRegion.Ann.rda")

eQTLsRanges <- GRanges(eQTLs.PUTM$chr, IRanges(eQTLs.PUTM$start, eQTLs.PUTM$end))

# we get overlap for each region
over <- lapply(caudate, function(x){overlapsAny(eQTLsRanges,x)})


tmp <- cbind(ann=eQTLs.PUTM$ourAnn,
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})


# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann=eQTLs.PUTM$ourAnn,marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp)

sapply(seq_along(summOver), function(i){
  barplot(t(summOver[[i]][,2]),beside = T,col = 1:2,main=names(summOver[i]),
          legend.text = c("Non-overlapping","Overlapping"),args.legend = list(x = "topleft"),
          names.arg = summOver[[i]][,1])})


sapply(seq_along(summOver), function(i){
  barplot(summOver[[i]][,2][,2]/apply(summOver[[i]][,2],1,sum),beside = T,col = 1:4,main=names(summOver[i]),
          names.arg = summOver[[i]][,1],ylab = "ratio")})


tmp1 <- as.matrix(as.data.frame(lapply(summOver,function(x){x[,2][,"TRUE"]/apply(x[,2],1,sum)})))
rownames(tmp1) <- summOver[[1]][,1]



load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.PUTM")
rm(covs,PUTM,SNIG)
genes <- rownames(RPKM.cqn)
rm(RPKM.cqn)
exonDef <- read.csv("data/general/transcriptomeInfo.csv")
exonDef <- exonDef[which(as.character(genes) %in% as.character(exonDef$names)),]
rm(genes)

genes <- unlist(lapply(strsplit(as.character(exonDef$names),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))
geneDef <- read.delim("data/general/ensemblGenes.txt")

exonDef <- cbind(exonDef,Chromosome.Name=geneDef[match(as.character(genes),as.character(geneDef$Ensembl.Gene.ID)),"Chromosome.Name"])

exonRanges <- GRanges(paste0("chr",exonDef$Chromosome.Name), IRanges(exonDef$start, exonDef$end))

library("biomaRt")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","gene_biotype"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=as.character(genes), mart=ensembl)

#head(geneNames)
exonDef <- cbind(exonDef,biotype=geneNames[match(as.character(genes),as.character(geneNames$ensembl_gene_id)),"gene_biotype"])
bp <- barplot(table(exonDef$biotype),main="Gene biotype for the brain expressed exons - PUTM")
text(bp , 0,labels = table(exonDef$biotype),cex=1,pos=3)

rm(genes,bp,geneDef)

rm(over)
over <- lapply(caudate, function(x){overlapsAny(exonRanges,x)})

tmp <- cbind(biotype=as.character(exonDef$biotype),
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})



## get the overlap based on the biotype

tmp <- cbind(tmp,allMarks=over$allMarks)

summOver <- list(DNAse=aggregate(DNase ~ biotype,data=tmp,table),
                 H3K27ac=aggregate(H3K27ac ~ biotype,data=tmp,table),
                 H3K4me1=aggregate(H3K4me1 ~ biotype,data=tmp,table),
                 H3K4me3=aggregate(H3K4me3 ~ biotype,data=tmp,table),
                 allMarks=aggregate(allMarks ~ biotype,data=tmp,table))

#' We now visualise the overlap for the brain annotated expressed exons divided by biotype 
barplot(sapply(seq_along(summOver), function(i){summOver[[i]][,2][,"TRUE"]/
    apply(summOver[[i]][,2],1,sum)}),
    names.arg = names(summOver),ylab = "ratio",
    main="Overlap annotated expressed exons by biotype - PUTM",beside = T,
    legend.text = summOver[[1]][,1],col=1:6,ylim=c(0,0.7),args.legend = list(x = "topleft"))


rm(summOver)

summOver <- lapply(over,table)


barplot(sapply(seq_along(summOver), function(i){summOver[[i]]["TRUE"]/sum(summOver[[i]])}),
        names.arg = names(summOver),ylab = "ratio", main="Overlap with expressed exons PUTM")

tmp1 <- as.matrix(rbind(tmp1,annExons=as.data.frame(lapply(summOver, function(x){x["TRUE"]/sum(x)}))[1,]))

barplot(tmp1,beside = T, 
        args.legend = list(x = "topleft"),legend.text = rownames(tmp1),col=1:5,ylab = "ratio",main="Overlap PUTM")



# library(GenomicFeatures)
# 
# write.table(as.data.frame(eQTLsRanges)[,1:3],col.names = F ,row.names = F,file = "tmp/eQTLs.SNIG.BED",quote = F)


## We check the overlap for SNPs and histone marks
{
#' we now check whether there is enrichement. 
rm(list=ls())
suppressMessages(library(devtools))
suppressMessages(library(GenomicRanges))
setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()
library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
opts_knit$set(root.dir=dirname(path$networkPathname))
rm(path)


H3K4me1 <- read.delim("data/Roadmap/E074-H3K4me1.imputed.gappedPeak.bed.gPk",header=F)
H3K4me3 <- read.delim("data/Roadmap/E074-H3K4me3.imputed.gappedPeak.bed.gPk",header=F)
H3K27ac <- read.delim("data/Roadmap/E074-H3K27ac.imputed.gappedPeak.bed.gPk",header=F)
DNase <- read.delim("data/Roadmap/E074-DNase.imputed.gappedPeak.bed.gPk",header=F)

#' convert all the ranges in Genomics range object
DNase <- GRanges(DNase[,1], IRanges(DNase[,2], DNase[,3]))
H3K4me1 <- GRanges(H3K4me1[,1], IRanges(H3K4me1[,2], H3K4me1[,3]))
H3K4me3 <- GRanges(H3K4me3[,1], IRanges(H3K4me3[,2], H3K4me3[,3]))
H3K27ac <- GRanges(H3K27ac[,1], IRanges(H3K27ac[,2], H3K27ac[,3]))

nigra <- list(DNase=DNase,H3K27ac=H3K27ac,H3K4me3=H3K4me3,H3K4me1=H3K4me1)
rm(H3K27ac,H3K4me3,H3K4me1,DNase)

#' So first thing I checked how much was the overlap, for example between two diffent histone marks for methylation. 
load("data/results/novelIntergenicRegionsFinal/intergenicRegion.Ann.rda")

chr <- unlist(lapply(strsplit(as.character(eQTLs.SNIG$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLs.SNIG$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann=eQTLs.SNIG$ourAnn,
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})

# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann=eQTLs.SNIG$ourAnn,marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp,pos,chr)



barplot(sapply(seq_along(summOver), function(i){summOver[[i]][,2][,"TRUE"]/
    apply(summOver[[i]][,2],1,sum)}),
    names.arg = names(summOver),ylab = "ratio",
    main="Overlap histones and eQTL Snps - SNIG",beside = T,
    legend.text = summOver[[1]][,1],col=1:4,ylim=c(0,0.7),args.legend = list(x = "topleft"))
head(summOver)

#' same thing for caudate
rm(list=ls())

H3K4me1 <- read.delim("data/Roadmap/E068-H3K4me1.imputed.gappedPeak.bed.gPk",header=F)
H3K4me3 <- read.delim("data/Roadmap/E068-H3K4me3.imputed.gappedPeak.bed.gPk",header=F)
H3K27ac <- read.delim("data/Roadmap/E068-H3K27ac.imputed.gappedPeak.bed.gPk",header=F)
DNase <- read.delim("data/Roadmap/E068-DNase.imputed.gappedPeak.bed.gPk",header=F)


#' convert all the ranges in Genomics range object
DNase <- GRanges(DNase[,1], IRanges(DNase[,2], DNase[,3]))
H3K4me1 <- GRanges(H3K4me1[,1], IRanges(H3K4me1[,2], H3K4me1[,3]))
H3K4me3 <- GRanges(H3K4me3[,1], IRanges(H3K4me3[,2], H3K4me3[,3]))
H3K27ac <- GRanges(H3K27ac[,1], IRanges(H3K27ac[,2], H3K27ac[,3]))

caudate <- list(DNase=DNase,H3K27ac=H3K27ac,H3K4me3=H3K4me3,H3K4me1=H3K4me1)
rm(H3K27ac,H3K4me3,H3K4me1,DNase)

load("data/results/novelIntergenicRegionsFinal/intergenicRegion.Ann.rda")

chr <- unlist(lapply(strsplit(as.character(eQTLs.PUTM$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLs.PUTM$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

# we get overlap for each region
over <- lapply(caudate, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann=eQTLs.PUTM$ourAnn,
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})

# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann=eQTLs.PUTM$ourAnn,marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp)



barplot(sapply(seq_along(summOver), function(i){summOver[[i]][,2][,"TRUE"]/
    apply(summOver[[i]][,2],1,sum)}),
    names.arg = names(summOver),ylab = "ratio",
    main="Overlap histones and eQTL Snps - PUTM",beside = T,
    legend.text = summOver[[1]][,1],col=1:4,ylim=c(0,0.7),args.legend = list(x = "topleft"))
head(summOver)



### we look for overlap of histone modification and annotated exons
H3K4me1 <- read.delim("data/Roadmap/E074-H3K4me1.imputed.gappedPeak.bed.gPk",header=F)
H3K4me3 <- read.delim("data/Roadmap/E074-H3K4me3.imputed.gappedPeak.bed.gPk",header=F)
H3K27ac <- read.delim("data/Roadmap/E074-H3K27ac.imputed.gappedPeak.bed.gPk",header=F)
DNase <- read.delim("data/Roadmap/E074-DNase.imputed.gappedPeak.bed.gPk",header=F)

#' convert all the ranges in Genomics range object
DNase <- GRanges(DNase[,1], IRanges(DNase[,2], DNase[,3]))
H3K4me1 <- GRanges(H3K4me1[,1], IRanges(H3K4me1[,2], H3K4me1[,3]))
H3K4me3 <- GRanges(H3K4me3[,1], IRanges(H3K4me3[,2], H3K4me3[,3]))
H3K27ac <- GRanges(H3K27ac[,1], IRanges(H3K27ac[,2], H3K27ac[,3]))

nigra <- list(DNase=DNase,H3K27ac=H3K27ac,H3K4me3=H3K4me3,H3K4me1=H3K4me1)
rm(H3K27ac,H3K4me3,H3K4me1,DNase)

#' So first thing I checked how much was the overlap, for example between two diffent histone marks for methylation. 
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")

chr <- unlist(lapply(strsplit(as.character(eQTLSNIG$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLSNIG$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann="geneExonic",
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})



# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann="geneExonic",marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp,pos,chr)

summOverAll <- summOver

#' So first thing I checked how much was the overlap, for example between two diffent histone marks for methylation. 
eQTLSNIG <- read.delim("data/results/finaleQTLs/exons.SNIG.txt",sep=" ")

chr <- unlist(lapply(strsplit(as.character(eQTLSNIG$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLSNIG$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann="exons",
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})



# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann="exons",marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp,pos,chr)

for(i in names(summOverAll))
{
  summOverAll[[i]] <- rbind(summOverAll[[i]],summOver[[i]])
}


load("data/results/finaleQTLs/eQTL.ExExJun.SNIG.rda")
eQTLSNIG <- eQTL.SNIG
rm(eQTL.SNIG)

chr <- unlist(lapply(strsplit(as.character(eQTLSNIG$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLSNIG$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann="exonsJun",
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})



# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann="exonsJun",marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp,pos,chr)

for(i in names(summOverAll))
{
  summOverAll[[i]] <- rbind(summOverAll[[i]],summOver[[i]])
}



barplot(sapply(seq_along(summOverAll), function(i){summOverAll[[i]][,2][,"TRUE"]/
    apply(summOverAll[[i]][,2],1,sum)}),
    names.arg = names(summOverAll),ylab = "ratio",
    main="Overlap histones and eQTL Snps - SNIG",beside = T,
    legend.text = summOverAll[[1]][,1],col=1:3,ylim=c(0,0.7),args.legend = list(x = "topleft"))


# We will do it for PUTM as well

rm(list=ls())

### we look for overlap of histone modification and annotated exons
H3K4me1 <- read.delim("data/Roadmap/E068-H3K4me1.imputed.gappedPeak.bed.gPk",header=F)
H3K4me3 <- read.delim("data/Roadmap/E068-H3K4me3.imputed.gappedPeak.bed.gPk",header=F)
H3K27ac <- read.delim("data/Roadmap/E068-H3K27ac.imputed.gappedPeak.bed.gPk",header=F)
DNase <- read.delim("data/Roadmap/E068-DNase.imputed.gappedPeak.bed.gPk",header=F)

#' convert all the ranges in Genomics range object
DNase <- GRanges(DNase[,1], IRanges(DNase[,2], DNase[,3]))
H3K4me1 <- GRanges(H3K4me1[,1], IRanges(H3K4me1[,2], H3K4me1[,3]))
H3K4me3 <- GRanges(H3K4me3[,1], IRanges(H3K4me3[,2], H3K4me3[,3]))
H3K27ac <- GRanges(H3K27ac[,1], IRanges(H3K27ac[,2], H3K27ac[,3]))

nigra <- list(DNase=DNase,H3K27ac=H3K27ac,H3K4me3=H3K4me3,H3K4me1=H3K4me1)
rm(H3K27ac,H3K4me3,H3K4me1,DNase)

#' So first thing I checked how much was the overlap, for example between two diffent histone marks for methylation. 
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")

chr <- unlist(lapply(strsplit(as.character(eQTLPUTM$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLPUTM$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann="geneExonic",
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})



# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann="geneExonic",marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp,pos,chr)

summOverAll <- summOver

#' So first thing I checked how much was the overlap, for example between two diffent histone marks for methylation. 
eQTLPUTM <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")

chr <- unlist(lapply(strsplit(as.character(eQTLPUTM$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLPUTM$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann="exons",
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})



# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann="exons",marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp,pos,chr)

for(i in names(summOverAll))
{
  summOverAll[[i]] <- rbind(summOverAll[[i]],summOver[[i]])
}


load("data/results/finaleQTLs/eQTL.ExExJun.PUTM.rda")
eQTLPUTM <- eQTL.PUTM
rm(eQTL.PUTM)

chr <- unlist(lapply(strsplit(as.character(eQTLPUTM$snps),":"),function(x){c(x[1])}))
pos <- unlist(lapply(strsplit(as.character(eQTLPUTM$snps),":"),function(x){c(x[2])}))

eQTLsRanges <- GRanges(chr, IRanges(as.numeric(pos), as.numeric(pos)))

over <- lapply(nigra, function(x){overlapsAny(eQTLsRanges,x)})

tmp <- cbind(ann="exonsJun",
             DNase=over$DNase,
             H3K27ac=over$H3K27ac,
             H3K4me3=over$H3K4me3,
             H3K4me1=over$H3K4me1)


#' We also get the overlap for all of them. We count whether at least the region overlap with at least one mark
over$allMarks <- apply(tmp,1,function(x){
  if(table(x=="TRUE")["FALSE"]<5){
    return(TRUE)  
  }else{
    return(FALSE)
  }
})



# we get the summary of the overlap
summOver <- lapply(over, function(x){
  tmp <- cbind(ann="exonsJun",marker=x)
  aggregate(marker ~ ann,data=tmp,table) 
})
rm(tmp,pos,chr)

for(i in names(summOverAll))
{
  summOverAll[[i]] <- rbind(summOverAll[[i]],summOver[[i]])
}



barplot(sapply(seq_along(summOverAll), function(i){summOverAll[[i]][,2][,"TRUE"]/
    apply(summOverAll[[i]][,2],1,sum)}),
    names.arg = names(summOverAll),ylab = "ratio",
    main="Overlap histones and eQTL Snps - PUTM",beside = T,
    legend.text = summOverAll[[1]][,1],col=1:3,ylim=c(0,0.7),args.legend = list(x = "topleft"))







}




