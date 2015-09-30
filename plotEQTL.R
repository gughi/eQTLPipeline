library("biomaRt")
library("Gviz")
## we select the archive version
listMarts(host="Jun2013.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

load("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr12.rda")
load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
#geneEQTL <- read.delim("data/results/genic/geneExons/fullResults/PUTM/ENSG00000204625")
geneEQTL <- read.delim("data/results/genic/geneIntronic/fullResults/PUTM/ENSG00000050426")

# eQTL in the in gene exons
chr <- 12
# start <- 29942889   
# end <- 29946183
gen = "hg19"
# leadingSNP <- "chr6:29955809"
     
# eQTL in the in gene intronic
# start <- 29942889   
# end <- 29946183
leadingSNP <- "chr12:52391833"




corLeadingSNP <- unlist(strsplit(leadingSNP,":"))[2]
start <- as.numeric(corLeadingSNP) -100000
end <- as.numeric(corLeadingSNP) +50000

fullCovtmp <- fullCov$chr12[(start):(end),PUTM$A.CEL_file]
dim(fullCovtmp)
rownames(fullCovtmp) <- (start):(end)



##tail(data.frame(fullCovtmp))

filters <- c("chromosomal_region")


chromoReg <- paste0(chr,":",start,":",end,":-1,",
                    chr,":",start,":",end,":1")


annReg <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=list(chromoReg), mart=ensembl)


width <- (annReg$exon_chrom_end - annReg$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
annReg <- cbind(annReg[,1:3],width,annReg[,4:9])  
rm(width)

annReg$strand <- gsub("-1","-",annReg$strand)
annReg$strand <- gsub("1","+",annReg$strand)
annReg$strand <- as.factor(annReg$strand)
## we now change the column names

colnames(annReg) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  



grtrack <- GeneRegionTrack(annReg, genome = gen,
                           chromosome = chr)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

## we load the dat
meanCov <- apply(data.frame(fullCovtmp),1,mean)
dat <- data.frame(paste0('chr',chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",type="histogram",name="coverage")
dtrack_heat <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",type="heatmap")

coords <- unlist(lapply(strsplit(as.character(geneEQTL$SNP),":"),function(x){x[2]}))


pvals <- as.data.frame(-log(geneEQTL$p.value))
pvals <- cbind(pvals,pvals)
colnames(pvals)<- c("red","black")

pvals[-match(corLeadingSNP,coords),1] <- NA
pvals[match(corLeadingSNP,coords),2] <- NA

pvalstrack <- DataTrack(data=pvals, start = as.numeric(coords),
                        end = (as.numeric(coords)+1), chromosome = chr, genome = gen,
                        name = "Pvalues",type=c("p"),groups=c("red","black"),col=c("black","red"),cex=1)


plotTracks(list(dtrack_heat,gtrack,dtrack,grtrack,pvalstrack), from = start,
           to = end, chromosome = chr, 
           transcriptAnnotation = "symbol")




## now I plot the eQTL that is in the intrnic and non in the exonic

## eQT
par(mfrow=c(1,2))

gene <- "ENSG00000204625"
snp <- "chr6:29915975:CA_C"
##snp <- "chr6:29955809"
dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

#exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/guelfi/plinkOutput/eigenvec"
title <- paste0("Gene Intronic ","(",snp,")")

snp <- gsub("chr", "", snp)
snps  <- unlist( read.table.rows(paste0(dosageFile), keepRows=snp, sep=" ") )
info  <- read.table.rows(infoFile, keepRows=snp, sep=" ", colClasses="character")


## load the expresssion
load(paste0(exprFile))

IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
PUTMtmp <- resids[,as.character(gene)]
rm(IDs)
PUTMsnp <- snps[names(PUTMtmp)]

identical(names(PUTMtmp),names(PUTMsnp))

## we get the genetic PCAs
my.covTMP <- read.table.rows(paste0(eigenFile), keepRows=names(PUTMtmp), sep=" ",header=F)
my.cov0 <- as.matrix(my.covTMP[names(PUTMsnp),2:4])
my.cov0 <- t(my.cov0)
rm(my.covTMP)

identical(names(PUTMtmp),colnames(my.cov0))

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="PUTM", main=title)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval",fit["PUTMsnp", "Pr(>|t|)"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )














