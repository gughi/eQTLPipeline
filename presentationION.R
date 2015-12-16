## For presentation of IoN
load("data/results/finaleQTLs/eQTL.ExExJun.PUTM.rda")

eQTL.PUTM[which(eQTL.PUTM$geneID %in% "ENSG00000197386"),]

gene <- "ENSG00000197386"
snp <- "chr4:2302211"
##

load(paste0("/home/seb/eQTL/snps/byGene/",gene,".rda"))
markers <- markers[snp,]


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

library(devtools)
load_all()

plotLoceQTLs(gene = gene,ensembl = ensembl,IDs = IDs,genotype = genotype,
             highLight=eQTL.PUTM[intersect(which(eQTL.PUTM$snps %in% snp),which(eQTL.PUTM$geneID %in% gene)),
                                 c("chrExon1","startExon1","endExon1","startExon2","endExon2","exExID")])



load(file="data/expr/normalisedCounts/genicIntergenic.rda")


RPKM.cqn <- RPKM.cqn[grep("DER",rownames(RPKM.cqn)),]


load("data/general/sampleInfo.rda")
# PCA PC1 vs PC2 

PCAres<- prcomp(t(RPKM.cqn),scale.=TRUE)
par(mfrow=c(1,1))
PUTM <- sampleInfo[sampleInfo[, 6] == "PUTM",]
SNIG <- sampleInfo[sampleInfo[, 6] == "SNIG",]

plot(PCAres, main="PCA axes intergenic quantification (PUTM + SNIG)")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue (intergenic)",ylab="",xlab="",ylim=c(-350,250), xlim=c(-450,850))
points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomright", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    

rm(PUTM,SNIG,PCAres)


load("data/expr/rawCounts/intergenic/allChromosomes.rda")

PCA <- prcomp(coverage)


rm(list=ls())

load("data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
eQTLPUTMEI <- eQTLPUTM
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTMI <- eQTLPUTM
load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")
eQTLPUTMInter <- eQTLPUTM
eQTLPUTMExons <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTMEx <- eQTLPUTM
load("data/results/finaleQTLs/eQTL.ExExJun.PUTM.rda")
eQTLPUTMExExJun <- eQTL.PUTM
rm(eQTLPUTM,eQTL.PUTM)
par(mar=c(3,2,3,1))

numeQTLs10 <- c(
                ##length(unique(eQTLPUTMEI$gene)),
                length(unique(eQTLPUTMEx$gene)),
                ##length(unique(eQTLPUTMI$gene)),
                length(unique(eQTLPUTMExons$gene)),
                length(unique(eQTLPUTMExExJun$exExID)),
                length(unique(eQTLPUTMInter$gene)))

names(numeQTLs10) <- c("Gene",
                      ##"GeneIntronic",
                      "Exons","exonExonJunc","Intergenic")

barplot(numeQTLs10,col='red',main= "eQTLs PUTM",border=F,
        sub=paste("PUTM"))

numeQTLs5 <- c(
              ##length(unique(eQTLPUTMEI[which(eQTLPUTMEI$myFDR<0.05),"gene"])),
               length(unique(eQTLPUTMEx[which(eQTLPUTMEx$myFDR<0.05),"gene"])),
               ##length(unique(eQTLPUTMI[which(eQTLPUTMI$myFDR<0.05),"gene"])),
               length(unique(eQTLPUTMExons[which(eQTLPUTMExons$myFDR<0.05),"gene"])),
               length(unique(eQTLPUTMExExJun[which(eQTLPUTMExExJun$myFDR<0.05),"exExID"])),
               length(unique(eQTLPUTMInter[which(eQTLPUTMInter$myFDR<0.05),"gene"])))

names(numeQTLs5) <- c("Gene",
                      ##"GeneIntronic",
                      "Exons","exonExonJunc","Intergenic")
barplot(numeQTLs5,add=T,col='orange',border=F)

numeQTLs1 <- c(
              ##length(unique(eQTLPUTMEI[which(eQTLPUTMEI$myFDR<0.01),"gene"])),
               length(unique(eQTLPUTMEx[which(eQTLPUTMEx$myFDR<0.01),"gene"])),
               ##length(unique(eQTLPUTMI[which(eQTLPUTMI$myFDR<0.01),"gene"])),
               length(unique(eQTLPUTMExons[which(eQTLPUTMExons$myFDR<0.01),"gene"])),
               length(unique(eQTLPUTMExExJun[which(eQTLPUTMExExJun$myFDR<0.01),"exExID"])),
               length(unique(eQTLPUTMInter[which(eQTLPUTMInter$myFDR<0.01),"gene"])))

names(numeQTLs1) <- c("Gene",
                      ##"GeneIntronic",
                      "Exons","exonExonJunc","Intergenic")

barplot(numeQTLs1,add=T,col='darkgreen',border=F)
legend("topleft",c("10%","5%","1%"),col=c('red','orange','darkgreen'),title="FDR",pch=15)

## table of the independecy of the signals
sign <- cbind(table(eQTLPUTM$degree),table(eQTLSNIG$degree))

par(mar=c(9,3,1,1))
barplot(sort(table(eQTLPUTM$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLPUTM$gene_biotype))),las=2)

barplot(sort(table(eQTLSNIG$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLSNIG$gene_biotype))),las=2,main="test")

rm(numeQTLs1,numeQTLs10,numeQTLs5)



rm(list=ls())


load("data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
eQTLSNIGEI <- eQTLSNIG
load("data/results/finaleQTLs/intronic.Ann.SNIG.rda")
eQTLSNIGI <- eQTLSNIG
load("data/results/finaleQTLs/intergenic.Ann.SNIG.rda")
eQTLSNIGInter <- eQTLSNIG
eQTLSNIGExons <- read.delim("data/results/finaleQTLs/exons.SNIG.txt",sep=" ")
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
eQTLSNIGEx <- eQTLSNIG
load("data/results/finaleQTLs/eQTL.ExExJun.SNIG.rda")
eQTLSNIGExExJun <- eQTL.SNIG
rm(eQTLSNIG,eQTL.SNIG)
par(mar=c(3,2,3,1))

numeQTLs10 <- c(
                ##length(unique(eQTLSNIGEI$gene)),
                length(unique(eQTLSNIGEx$gene)),
                ##length(unique(eQTLSNIGI$gene)),
                length(unique(eQTLSNIGExons$gene)),
                length(unique(eQTLSNIGExExJun$exExID)),
                length(unique(eQTLSNIGInter$gene)))

names(numeQTLs10) <- c("Gene",
                      ##"GeneIntronic",
                      "Exons","exonExonJunc","Intergenic")

barplot(numeQTLs10,col='red',main= "eQTLs SNIG",border=F,
        sub=paste("SNIG"))

numeQTLs5 <- c(
                ##length(unique(eQTLSNIGEI[which(eQTLSNIGEI$myFDR<0.05),"gene"])),
               length(unique(eQTLSNIGEx[which(eQTLSNIGEx$myFDR<0.05),"gene"])),
               ##length(unique(eQTLSNIGI[which(eQTLSNIGI$myFDR<0.05),"gene"])),
               length(unique(eQTLSNIGExons[which(eQTLSNIGExons$myFDR<0.05),"gene"])),
               length(unique(eQTLSNIGExExJun[which(eQTLSNIGExExJun$myFDR<0.05),"exExID"])),
               length(unique(eQTLSNIGInter[which(eQTLSNIGInter$myFDR<0.05),"gene"])))

names(numeQTLs5) <- c("Gene",
                      ##"GeneIntronic",
                      "Exons","exonExonJunc","Intergenic")

barplot(numeQTLs5,add=T,col='orange',border=F)

numeQTLs1 <- c(##length(unique(eQTLSNIGEI[which(eQTLSNIGEI$myFDR<0.01),"gene"])),
               length(unique(eQTLSNIGEx[which(eQTLSNIGEx$myFDR<0.01),"gene"])),
               ##length(unique(eQTLSNIGI[which(eQTLSNIGI$myFDR<0.01),"gene"])),
               length(unique(eQTLSNIGExons[which(eQTLSNIGExons$myFDR<0.01),"gene"])),
               length(unique(eQTLSNIGExExJun[which(eQTLSNIGExExJun$myFDR<0.01),"exExID"])),
               length(unique(eQTLSNIGInter[which(eQTLSNIGInter$myFDR<0.01),"gene"])))

names(numeQTLs1) <- c("Gene",
                      ##"GeneIntronic",
                      "Exons","exonExonJunc","Intergenic")

barplot(numeQTLs1,add=T,col='darkgreen',border=F)
legend("topleft",c("10%","5%","1%"),col=c('red','orange','darkgreen'),title="FDR",pch=15)



head(eQTLSNIGEx)

###############
### TSS #######
###############

library("biomaRt")

eQTLPUTM <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

genes <- unlist(lapply(strsplit(as.character(eQTLPUTM$gene),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=genes, mart=ensembl)


load_all()
TSS <- sapply(genes, function(x){getTSS(x,geneNames)})

eQTLPUTM <- cbind(eQTLPUTM,TSS)

posSNP <- unlist(lapply(strsplit(as.character(eQTLPUTM$snps),":"),function(x){x[2]}))
DisGeneStart <- eQTLPUTM$TSS - as.integer(posSNP)
eQTLPUTM <- cbind(eQTLPUTM,DisGeneStart) 

eQTLPUTM[which(eQTLPUTM$DisGeneStart %in% min(eQTLPUTM$DisGeneStart)),]
hist(eQTLPUTM$DisGeneStart,main="Distance from TSS exon eQTLs PUTM",breaks=30)

max(table(genes)[which((table(genes)>1) %in% TRUE)])

genTmp <- table(genes)[which((table(genes)>1) %in% TRUE)] 
hist(eQTLPUTM[which(genes %in% names(genTmp)),"DisGeneStart"],breaks=40,probability=T)


library("biomaRt")





eQTLSNIG <- read.delim("data/results/finaleQTLs/exons.SNIG.txt",sep=" ")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

genes <- unlist(lapply(strsplit(as.character(eQTLSNIG$gene),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=genes, mart=ensembl)

TSS <- sapply(genes, function(x){getTSS(x,geneNames)})

eQTLSNIG <- cbind(eQTLSNIG,TSS)

posSNP <- unlist(lapply(strsplit(as.character(eQTLSNIG$snps),":"),function(x){x[2]}))
DisGeneStart <- eQTLSNIG$TSS - as.integer(posSNP)
eQTLSNIG <- cbind(eQTLSNIG,DisGeneStart) 

eQTLSNIG[which(eQTLSNIG$DisGeneStart %in% min(eQTLSNIG$DisGeneStart)),]
hist(eQTLSNIG$DisGeneStart,main="Distance from TSS exon eQTLs SNIG",breaks=30)

genTmp <- table(genes)[which((table(genes)>1) %in% TRUE)] 
hist(eQTLSNIG[which(genes %in% names(genTmp)),"DisGeneStart"],breaks=30)



rm(list=ls())

head(conse)

load("data/results/finaleQTLs/geneExonic.un.SNPAnn.PUTM.rda")


tmp <- sapply(conse[,2],function(x){unlist(strsplit(x,";"))[1]})

par(mar=c(11,2,3,1))
barplot(sort(table(tmp),decreasing=T),las=2, 
        main="Variants annotation PUTM",
        yaxt='n',col=2:length(table(tmp)))



load("data/results/finaleQTLs/geneExonic.un.rda")
tmp <- sapply(conse[,2],function(x){unlist(strsplit(x,";"))[1]})

par(mar=c(11,2,3,1))
barplot(sort(table(tmp),decreasing=T),las=2, 
        main="Variants annotation SNIG",
        yaxt='n',col=2:length(table(tmp)))



load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")
load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")

tmpExpr <- apply(RPKM.cqn,1,mean)

tmpExpr <- tmpExpr[as.character(eQTLPUTM$gene)]
starStopReg <- starStopReg[names(sort(tmpExpr,decreasing=T)),]

tmpExpr <- sort(tmpExpr,decreasing=T)

ensemblDef <- read.delim("data/general/ensemblGenes.txt",row.names=1)

tmppp <- GRanges(paste0("chr",ensemblDef[,3]), IRanges(ensemblDef[,1], ensemblDef[,2]))
names(tmppp) <- rownames(ensemblDef)
ensemblDef <- tmppp
rm(tmppp)
genDef <- GRanges(paste0(starStopReg[,1]), IRanges(starStopReg[,2], starStopReg[,3]))
names(genDef) <- rownames(starStopReg)

head(tmpExpr[which(countOverlaps(genDef, ensemblDef)==0)],20)

(genDef["DER2044",], ensemblDef)


derReg <- "DER2044"

RPKM.cqn["DER1885",]
# eQTL in the in gene exons
chr <- gsub("chr","",starStopReg[derReg,1])
start <- starStopReg[derReg,2] -10000
end <- starStopReg[derReg,3]   +10000
gen = "hg19"
leadingSNP <- eQTLPUTM[which(eQTLPUTM$gene %in% derReg),1]


corLeadingSNP <- unlist(strsplit(as.character(leadingSNP),":"))[2]
#start <- as.numeric(corLeadingSNP) -10000
#end <- as.numeric(corLeadingSNP) +50000


load(paste0("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr",chr,".rda"))
load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
geneEQTL <- read.delim(paste0("data/results/intergenic/fullResults/PUTM/",derReg))


fullCovtmp <- fullCov[[1]][(start):(end),PUTM$A.CEL_file]
dim(fullCovtmp)
rownames(fullCovtmp) <- (start):(end)
##tail(data.frame(fullCovtmp))

filters <- c("chromosomal_region")


chromoReg <- paste0(chr,":",start,":",end,":-1,",
                    chr,":",start,":",end,":1")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

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
#itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

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
pvals <- as.data.frame(-log(geneEQTL$FDR))
pvals <- cbind(pvals,pvals)
colnames(pvals)<- c("red","black")

pvals[-match(corLeadingSNP,coords),1] <- NA
pvals[match(corLeadingSNP,coords),2] <- NA

pvalstrack <- DataTrack(data=pvals, start = as.numeric(coords),
                        end = (as.numeric(coords)+1), chromosome = chr, genome = gen,
                        name = "FDR",type=c("p"),groups=c("red","black"),col=c("black","red"),cex=1)


plotTracks(list(gtrack,dtrack_heat,dtrack,grtrack,pvalstrack), from = start,
           to = end, chromosome = chr, 
           transcriptAnnotation = "symbol")












load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")

derReg <- "DER2044"

# eQTL in the in gene exons
chr <- gsub("chr","",starStopReg[derReg,1])
start <- starStopReg[derReg,2] -10000
end <- starStopReg[derReg,3]   +10000
gen = "hg19"



chromoReg <- paste0(chr,":",start,":",end,":-1,",
                    chr,":",start,":",end,":1")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")
filters <- c("chromosomal_region")

defGen <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=list(chromoReg), mart=ensembl)
rm(geneType,sta,filters)
width <- (defGen$exon_chrom_end - defGen$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
defGen <- cbind(defGen[,1:3],width,defGen[,4:9])  
rm(width)

defGen$strand <- gsub("-1","-",defGen$strand)
defGen$strand <- gsub("1","+",defGen$strand)
defGen$strand <- as.factor(defGen$strand)
## we now change the column names

colnames(defGen) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  


load(paste0("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr",chr,".rda"))

## we load the data for the example

### Select of the bp for the gene considered
load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]


fullCovtmp <- fullCov[[1]][(start):(end),as.character(PUTM$A.CEL_file)]
rownames(fullCovtmp) <- (start):(end)

tail(data.frame(fullCovtmp))


grtrack <- GeneRegionTrack(defGen, genome = gen,
                           chromosome = chr, name = "intergenic Region")

gtrack <- GenomeAxisTrack()
##itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

## we load the dat
snp <- eQTLPUTM[which(eQTLPUTM$gene %in% derReg),1]

load(paste0("data/snps/byRegion/PUTM/",derReg,".rda"))
markers <- markers[as.character(snp),]
IDs <- gsub("/","_",PUTM$U.SD_No)
tmp <- markers[,as.character(IDs)]
names(tmp) <- as.character(PUTM$A.CEL_file)
markers <- list(info=markers[,c(1:6)],genotype=tmp)




tmp <- round(as.numeric(markers$genotype))
## we load 3expression based on teh genotype
IDsGen <- names(markers$genotype[which(tmp %in% "0")]) 
meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
meanAll <- meanCov
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_gRef <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
rm(dat)


## we load 3expression based on teh genotype
IDsGen <- names(markers$genotype[which(tmp %in% "1")]) 
meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
meanAll <- cbind(meanAll,meanCov)
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_gHet <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
rm(dat)

## we load 3expression based on teh genotype
IDsGen <- names(markers$genotype[which(tmp %in% "2")]) 
rm(tmp)
meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
meanAll <- cbind(meanAll,meanCov)
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_gAlt <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
rm(dat)


colnames(meanAll) <- c(paste0(markers$info$Al1,markers$info$Al1),
                       paste0(markers$info$Al1,markers$info$Al2),
                       paste0(markers$info$Al2,markers$info$Al2))


library(MatrixEQTL)
my.expr <- SlicedData$new()
my.expr$CreateFromMatrix(as.matrix(as.data.frame(fullCovtmp[,names(markers$genotype)])))
my.markers <- SlicedData$new()
my.markers$CreateFromMatrix(as.matrix(markers$genotype))
store <- Matrix_eQTL_main( my.markers, my.expr, output_file_name = NULL,pvOutputThreshold=1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
head(store$all$eqtls)

## betas
datBetas <- GRanges(chr, IRanges(as.numeric(as.character(store$all$eqtls$gene))
                                 , as.numeric(as.character(store$all$eqtls$gene))), cov=store$all$eqtls$beta)

dtrackBetas <- DataTrack(range=datBetas,chromosome=paste0("chr",chr),genome="hg19",name="betas",type="histogram")

datPval <- GRanges(chr, IRanges(as.numeric(as.character(store$all$eqtls$gene))
                                , as.numeric(as.character(store$all$eqtls$gene))), cov=-log10(store$all$eqtls$pvalue))

dtrackPval <- DataTrack(range=datPval,chromosome=paste0("chr",chr),genome="hg19",name="-log10(pvalues)",type="gradient")                          


## The code below needs to improve
colnames(meanAll) <- c(paste0(markers$info$Al1,markers$info$Al1),paste0(markers$info$Al1,markers$info$Al2),paste0(markers$info$Al2,markers$info$Al2))
allInSameTrack <- DataTrack(data=t(meanAll),start=as.numeric(rownames(meanAll)),end=as.numeric(rownames(meanAll)), chromosome = chr, genome = gen,
                            name = "Stratified raw counts",type=c("l"),groups=c(paste0(markers$info$Al1,markers$info$Al1),
                                                                                paste0(markers$info$Al1,markers$info$Al2),
                                                                                paste0(markers$info$Al2,markers$info$Al2)),legend = TRUE
                            ,col=c("black","red","blue"))
#
  plotTracks(list(gtrack,dtrackPval,dtrackBetas,allInSameTrack),transcriptAnnotation = "symbol")



library (plyr)

gwasCat <- read.delim("data/general/GWAS_Disease_Classification_070415.tsv")
gwasCat <- gwasCat[rownames((unique(gwasCat[,1:9]))),]
my.eQTLs   <- read.delim(file="data/results/finaleQTLs/exons.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(my.eQTLs$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
my.eQTLs <- df[which(df[,7] < 0.05),]
dim(my.eQTLs)
hits <- gsub("chr","",my.eQTLs[,1])

cat(hits, sep="\n", file=(tmpf <- tempfile()))
map <- read.delim( pipe(paste0("fgrep -w --file=", tmpf, " /home/adai/rawData/chrpos2rsid_UKBEC")), sep=" ", header=FALSE, as.is=TRUE )
hits <- map$V2[!is.na(map$V2)]

hitsGwas <- unlist(lapply(strsplit(as.character(gwasCat$Strongest.SNP.Risk.Allele),"-"),function(x){x[1]}))

out <- gwasCat[which(hitsGwas %in% hits),]

print(nrow(out))
print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))


# > print(nrow(out))
# [1] 124
# > print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
# [1] "number of SNPs found in the GWAS catalogue related with brain: 12"
# > print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
# [1] "number of SNPs found in the GWAS catalogue related with Adult neurological disease: 8"
# > print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
# [1] "rate Adult Neuro Br: 0.0645161290322581"
# > print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))
# [1] "rate: 0.0967741935483871"


"chr17:46710944"
"chr7:100406823"

gwasCat[rownames((unique(gwasCat[,1:9]))),]

table(uniqueGWAs[,"General.Phenotype"]/)
978+7860

978/8838





library (plyr)

gwasCat <- read.delim("data/general/GWAS_Disease_Classification_070415.tsv")
gwasCat <- gwasCat[rownames((unique(gwasCat[,1:9]))),]
my.eQTLs   <- read.delim(file="data/results/finaleQTLs/exons.unsentinalised.PUTM.txt",  as.is=T, header=T)
tmp <- strsplit(my.eQTLs$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
my.eQTLs <- df[which(df[,7] < 0.05),]
dim(my.eQTLs)
hits <- gsub("chr","",my.eQTLs[,1])

cat(hits, sep="\n", file=(tmpf <- tempfile()))
map <- read.delim( pipe(paste0("fgrep -w --file=", tmpf, " /home/adai/rawData/chrpos2rsid_UKBEC")), sep=" ", header=FALSE, as.is=TRUE )
hits <- map$V2[!is.na(map$V2)]

hitsGwas <- unlist(lapply(strsplit(as.character(gwasCat$Strongest.SNP.Risk.Allele),"-"),function(x){x[1]}))

out <- gwasCat[which(hitsGwas %in% hits),]

print(nrow(out))
print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))

# > print(nrow(out))
# [1] 252
# > print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
# [1] "number of SNPs found in the GWAS catalogue related with brain: 45"
# > print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
# [1] "number of SNPs found in the GWAS catalogue related with Adult neurological disease: 15"
# > print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
# [1] "rate Adult Neuro Br: 0.0595238095238095"
# > print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))
# [1] "rate: 0.178571428571429"


library (plyr)

gwasCat <- read.delim("data/general/GWAS_Disease_Classification_070415.tsv")
gwasCat <- gwasCat[rownames((unique(gwasCat[,1:9]))),]
my.eQTLs   <- read.delim(file="data/results/finaleQTLs/eQTLgeneExons.unsentinalised.PUTM.txt",  as.is=T, header=T)
tmp <- strsplit(my.eQTLs$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
my.eQTLs <- df[which(df[,7] < 0.05),]
dim(my.eQTLs)
hits <- gsub("chr","",my.eQTLs[,1])

cat(hits, sep="\n", file=(tmpf <- tempfile()))
map <- read.delim( pipe(paste0("fgrep -w --file=", tmpf, " /home/adai/rawData/chrpos2rsid_UKBEC")), sep=" ", header=FALSE, as.is=TRUE )
hits <- map$V2[!is.na(map$V2)]

hitsGwas <- unlist(lapply(strsplit(as.character(gwasCat$Strongest.SNP.Risk.Allele),"-"),function(x){x[1]}))

out <- gwasCat[which(hitsGwas %in% hits),]

print(nrow(out))
print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))

# > print(nrow(out))
# [1] 72
# > print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
# [1] "number of SNPs found in the GWAS catalogue related with brain: 11"
# > print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
# [1] "number of SNPs found in the GWAS catalogue related with Adult neurological disease: 7"
# > print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
# [1] "rate Adult Neuro Br: 0.0972222222222222"
# > print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))
# [1] "rate: 0.152777777777778"
out[which(out$General.Phenotype %in% "brain"),]

library (plyr)

gwasCat <- read.delim("data/general/GWAS_Disease_Classification_070415.tsv")
gwasCat <- gwasCat[rownames((unique(gwasCat[,1:9]))),]
my.eQTLs   <- read.delim(file="data/results/finaleQTLs/eQTLgeneExons.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(my.eQTLs$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
my.eQTLs <- df[which(df[,7] < 0.05),]
dim(my.eQTLs)
hits <- gsub("chr","",my.eQTLs[,1])

cat(hits, sep="\n", file=(tmpf <- tempfile()))
map <- read.delim( pipe(paste0("fgrep -w --file=", tmpf, " /home/adai/rawData/chrpos2rsid_UKBEC")), sep=" ", header=FALSE, as.is=TRUE )
hits <- map$V2[!is.na(map$V2)]

hitsGwas <- unlist(lapply(strsplit(as.character(gwasCat$Strongest.SNP.Risk.Allele),"-"),function(x){x[1]}))

out <- gwasCat[which(hitsGwas %in% hits),]

print(nrow(out))
print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))

out[which(out$General.Phenotype %in% "brain"),c("Disease.Trait","ens72.chrpos")]
out

## "chr17:44856641"
## chr16:31095171 ENSG00000149923 PPP4C
## chr17:44788310 ENSG00000214401 
17:44828931
17:43513441

17:44788310     
17:43513441
7:100004446
my.eQTLs[which(my.eQTLs[,1]%in%"chr17:44788310"),]


# [,1]             [,2]                   [,3]                [,4]                   [,5]                   [,6]                 [,7]                  
# [1,] "chr17:44788310" "ENSG00000186868:E023" "-5.7039412623925"  "1.19422361056538e-07" "2.01079894332667e-07" "-0.421889923968665" "0.000674019805803099"
# [2,] "chr17:44788310" "ENSG00000214401:E006" "-5.28095404781991" "7.51368372363039e-07" "1.16352560844697e-06" "-0.710750049512539" "0.00395069490188486" 
# [3,] "chr17:44788310" "ENSG00000214401:E005" "-5.41467298016003" "4.23212259014192e-07" "6.57695614532827e-07" "-0.705055645668405" "0.00222525005789662" 
# [4,] "chr17:44788310" "ENSG00000186868:E022" "-5.86741564580864" "5.76486826647822e-08" "9.77965629576287e-08" "-0.211143672694628" "0.000325369164960031"
# [5,] "chr17:44788310" "ENSG00000120071:E006" "4.82440921831798"  "5.03968336869151e-06" "9.95865937688146e-06" "0.140498629862473"  "0.0287261952015416"  
# [6,] "chr17:44788310" "ENSG00000120071:E027" "5.51568773890372"  "2.73036494012197e-07" "5.10225254775129e-07" "0.232857762341719"  "0.00155630801586952" 
# [7,] "chr17:44788310" "ENSG00000120071:E036" "-4.76366521640029" "6.44643714718855e-06" "1.09294145565065e-05" "-0.392572484724611" "0.0367446917389747"  
# [8,] "chr17:44788310" "ENSG00000186868:E007" "-5.90228038540229" "4.92971879936354e-08" "8.35535522630865e-08" "-0.237739218118307" "0.000278233329036078"
# [9,] "chr17:44788310" "ENSG00000186868:E015" "-5.27294596125087" "7.77469989385629e-07" "1.31064534650313e-06" "-0.238230639541333" "0.00438804062009249" 
# [10,] "chr17:44788310" "ENSG00000186868:E009" "-5.64534529032552" "1.54693123668874e-07" "2.62019913136509e-07" "-0.238155382772076" "0.000873087989987122"
# 

transcriptomeInfo <- read.csv("data/general/transcriptomeInfo.csv")
head(transcriptomeInfo)
exonSel <- c("ENSG00000186868:E007","ENSG00000186868:E009","ENSG00000186868:E015","ENSG00000186868:E022","ENSG00000186868:E023" )
transcriptomeInfo[exonSel,]


gene <- "ENSG00000186868"

# eQTL in the in gene exons
chr <- "17"
start <- 44039704 - 1000
end <-  44067441+1000
gen = "hg19"

chromoReg <- paste0(chr,":",start,":",end,":-1,",
                    chr,":",start,":",end,":1")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")
filters <- c("chromosomal_region")

defGen <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=list(chromoReg), mart=ensembl)
rm(geneType,sta,filters)
width <- (defGen$exon_chrom_end - defGen$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
defGen <- cbind(defGen[,1:3],width,defGen[,4:9])  
rm(width)

defGen$strand <- gsub("-1","-",defGen$strand)
defGen$strand <- gsub("1","+",defGen$strand)
defGen$strand <- as.factor(defGen$strand)
## we now change the column names

colnames(defGen) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  


load(paste0("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr",chr,".rda"))

## we load the data for the example

### Select of the bp for the gene considered
load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]


fullCovtmp <- fullCov[[1]][(start):(end),as.character(PUTM$A.CEL_file)]
rownames(fullCovtmp) <- (start):(end)

tail(data.frame(fullCovtmp))


grtrack <- GeneRegionTrack(defGen, genome = gen,
                           chromosome = chr, name = "intergenic Region")

gtrack <- GenomeAxisTrack()
##itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

## we load the dat
snp <- "chr17:44788310"

load(paste0("/home/seb/eQTL/snps/byGene/",gene,".rda"))
markers <- markers[as.character(snp),]
IDs <- gsub("/","_",PUTM$U.SD_No)
tmp <- markers[,as.character(IDs)]
names(tmp) <- as.character(PUTM$A.CEL_file)
markers <- list(info=markers[,c(1:6)],genotype=tmp)




tmp <- round(as.numeric(markers$genotype))
## we load 3expression based on teh genotype
IDsGen <- names(markers$genotype[which(tmp %in% "0")]) 
meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
meanAll <- meanCov
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_gRef <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
rm(dat)


## we load 3expression based on teh genotype
IDsGen <- names(markers$genotype[which(tmp %in% "1")]) 
meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
meanAll <- cbind(meanAll,meanCov)
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_gHet <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
rm(dat)

## we load 3expression based on teh genotype
IDsGen <- names(markers$genotype[which(tmp %in% "2")]) 
rm(tmp)
meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
meanAll <- cbind(meanAll,meanCov)
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_gAlt <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
rm(dat)


colnames(meanAll) <- c(paste0(markers$info$Al1,markers$info$Al1),
                       paste0(markers$info$Al1,markers$info$Al2),
                       paste0(markers$info$Al2,markers$info$Al2))



library(MatrixEQTL)
my.expr <- SlicedData$new()
my.expr$CreateFromMatrix(as.matrix(as.data.frame(fullCovtmp[,names(markers$genotype)])))
my.markers <- SlicedData$new()
my.markers$CreateFromMatrix(as.matrix(markers$genotype))
store <- Matrix_eQTL_main( my.markers, my.expr, output_file_name = NULL,pvOutputThreshold=1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
head(store$all$eqtls)

## betas
datBetas <- GRanges(chr, IRanges(as.numeric(as.character(store$all$eqtls$gene))
                                 , as.numeric(as.character(store$all$eqtls$gene))), cov=store$all$eqtls$beta)

dtrackBetas <- DataTrack(range=datBetas,chromosome=paste0("chr",chr),genome="hg19",name="betas",type="histogram")

datPval <- GRanges(chr, IRanges(as.numeric(as.character(store$all$eqtls$gene))
                                , as.numeric(as.character(store$all$eqtls$gene))), cov=-log10(store$all$eqtls$pvalue))

dtrackPval <- DataTrack(range=datPval,chromosome=paste0("chr",chr),genome="hg19",name="-log10(pvalues)",type="gradient")                          


## The code below needs to improve
colnames(meanAll) <- c(paste0(markers$info$Al1,markers$info$Al1),paste0(markers$info$Al1,markers$info$Al2),paste0(markers$info$Al2,markers$info$Al2))
allInSameTrack <- DataTrack(data=t(meanAll),start=as.numeric(rownames(meanAll)),end=as.numeric(rownames(meanAll)), chromosome = chr, genome = gen,
                            name = "Stratified raw counts",type=c("l"),groups=c(paste0(markers$info$Al1,markers$info$Al1),
                                                                                paste0(markers$info$Al1,markers$info$Al2),
                                                                                paste0(markers$info$Al2,markers$info$Al2)),legend = TRUE
                            ,col=c("black","red","blue"))
#
plotTracks(list(gtrack,dtrackPval,dtrackBetas,allInSameTrack),transcriptAnnotation = "symbol",from = start, to = end)
plotTracks(list(gtrack,grtrack),transcriptAnnotation = "symbol")

plotTracks(list(gtrack,dtrackBetas,allInSameTrack),transcriptAnnotation = "symbol",from = start, to = end)

plotTracks(list(gtrack,dtrackPval,dtrackBetas,allInSameTrack),transcriptAnnotation = "symbol",from = 44067044, to = 44067541)




# ENSG00000186868:E007 44039704 44039836   133 ENSG00000186868:E007
# ENSG00000186868:E009 44049225 44049311    87 ENSG00000186868:E009
# ENSG00000186868:E015 44055741 44055746     6 ENSG00000186868:E015
# ENSG00000186868:E022 44064406 44064461    56 ENSG00000186868:E022
# ENSG00000186868:E023 44067244 44067441   198 ENSG00000186868:E023






eQTL.PUTM[which(eQTL.PUTM$geneID %in% "ENSG00000186868"),]

gene <- "ENSG00000186868"
snp <- "chr17:44788310"
##

load(paste0("/home/seb/eQTL/snps/byGene/",gene,".rda"))
markers <- markers[snp,]


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

library(devtools)
load_all()

plotLoceQTLs(gene = gene,ensembl = ensembl,IDs = IDs,genotype = genotype)




gene <- "ENSG00000186868:E007"
snp <- "chr17:44788310"
##snp <- "chr6:29955809"
dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

#exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
exprFile <- "data/expr/normalisedCounts/genic/exons/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0("MAPT","(",snp,")")

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






my.eQTLs <- read.delim(file="data/results/finaleQTLs/eQTLgeneExons.unsentinalised.PUTM.txt",  as.is=T, header=T)
tmp <- strsplit(my.eQTLs$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
my.eQTLs.PUTM <- df[which(df[,7] < 0.01),]
rm(my.eQTLs,df)

my.eQTLs <- read.delim(file="data/results/finaleQTLs/eQTLgeneExons.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(my.eQTLs$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
my.eQTLs.SNIG <- df[which(df[,7] < 0.01),]
rm(my.eQTLs,df)


PUTM <- paste0(my.eQTLs.PUTM[,1],my.eQTLs.PUTM[,2])
SNIG <- paste0(my.eQTLs.SNIG[,1],my.eQTLs.SNIG[,2])


length(intersect(PUTM,SNIG))/length(SNIG)


load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")


eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<=0.05),]
eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<=0.05),]
length(intersect(eQTLPUTM$gene,eQTLSNIG$gene))/length(unique(eQTLSNIG$gene))



load("data/results/finaleQTLs/eQTL.ExExJun.PUTM.rda")
load("data/results/finaleQTLs/eQTL.ExExJun.SNIG.rda")
eQTL.PUTM <- eQTL.PUTM[which(eQTL.PUTM$myFDR<=0.05),]
eQTL.SNIG <- eQTL.SNIG[which(eQTL.SNIG$myFDR<=0.05),]
length(intersect(eQTL.PUTM$geneID,eQTL.SNIG$geneID))/length(unique(eQTL.SNIG$geneID))



eQTL.PUTM <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")
eQTL.SNIG <- read.delim("data/results/finaleQTLs/exons.SNIG.txt",sep=" ")
eQTL.PUTM <- eQTL.PUTM[which(eQTL.PUTM$myFDR<=0.05),]
eQTL.SNIG <- eQTL.SNIG[which(eQTL.SNIG$myFDR<=0.05),]

genesPUTM <- unlist(lapply(strsplit(as.character(eQTL.PUTM$gene),":"),function(x){x[1]}))
genesPUTM <- unlist(sapply(strsplit(as.character(genesPUTM),"+",fixed=T),function(x){x[1]}))

genesSNIG <- unlist(lapply(strsplit(as.character(eQTL.SNIG$gene),":"),function(x){x[1]}))
genesSNIG <- unlist(sapply(strsplit(as.character(genesSNIG),"+",fixed=T),function(x){x[1]}))

length(intersect(genesPUTM,genesSNIG))/length(unique(genesSNIG))



chr <- 21
start <- 27543189
end <- 27589700
gen = "hg19"


filters <- c("chromosomal_region")
chromoReg <- paste0(chr,":",start,":",end,":-1,",
                    chr,":",start,":",end,":1")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

defGen <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=list(chromoReg), mart=ensembl)
chr <- unique(defGen$chromosome) 
rm(geneType,sta,filters)
width <- (defGen$exon_chrom_end - defGen$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
defGen <- cbind(defGen[,1:3],width,defGen[,4:9])  
rm(width)

defGen$strand <- gsub("-1","-",defGen$strand)
defGen$strand <- gsub("1","+",defGen$strand)
defGen$strand <- as.factor(defGen$strand)
## we now change the column names

colnames(defGen) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  


load(paste0("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr",chr,".rda"))

## we load the data for the example

### Select of the bp for the gene considered

if(!is.na(IDs))
{
  fullCovtmp <- fullCov[[1]][(startStop$start_position-10000):(startStop$end_position+10000),as.character(IDs)]
}else{
  fullCovtmp <- fullCov[[1]][(start-10000):(end+10000),]
}  
rownames(fullCovtmp) <- (start-10000):(end+10000)

tail(data.frame(fullCovtmp))

grtrack <- GeneRegionTrack(defGen, genome = gen,
                           chromosome = chr, name="Ensembl reference v72")
gtrack <- GenomeAxisTrack()
##itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

## we load the dat
meanCov <- apply(data.frame(fullCovtmp),1,mean)
dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",type="histogram",name="coverage")
dtrack_heat <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",type="heatmap")

## itrack <- IdeogramTrack(genome = gen, chromosome = chr)

plotTracks(list(gtrack, dtrack,grtrack),transcriptAnnotation = "symbol",from=27530000,to=27598000)                     












ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Feb2014.archive.ensembl.org",
                      dataset="hsapiens_gene_ensembl")
filters <- c("chromosomal_region")

chromoReg <- paste0(chr,":",start,":",27598000,":-1,",
                    chr,":",start,":",27598000,":1")


defGen <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=list(chromoReg), mart=ensembl)
chr <- unique(defGen$chromosome) 
rm(geneType,sta,filters)
width <- (defGen$exon_chrom_end - defGen$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
defGen <- cbind(defGen[,1:3],width,defGen[,4:9])  
rm(width)

defGen$strand <- gsub("-1","-",defGen$strand)
defGen$strand <- gsub("1","+",defGen$strand)
defGen$strand <- as.factor(defGen$strand)
## we now change the column names

colnames(defGen) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  

grtrack <- GeneRegionTrack(defGen, genome = gen,
                           chromosome = chr, name="Ensembl reference v75")

plotTracks(list(gtrack, dtrack,grtrack),transcriptAnnotation = "symbol",from=27530000,to=27598000)                     



