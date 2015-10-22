###  Manuel Sebastian Guelfi
### 21-5-15
### This script qunatify the abundance of the intronic region in each gene.

### the approach we will use is to load the expression quantification of the exonic regions only 
### produced by HTSeq-counts  and then we will load the expression of the genes accounting for
### exons and introns (quantified with intersect bed) and the we will subtract the first dataset.

rm(list=ls())
setwd("/home/guelfi/eQTLPipeline")
## load the reads for exonic regions
load("data/expr/rawCounts/genic/exprSQ.rda")
## dim(SQcounts)
## [1]   170 64078

## load the reads for intronic and exonic regions
SQEIcounts <- read.csv("data/expr/rawCounts/genic/exprExonIntr.csv",row.names=1)
## dim(SQEIcounts)
## [1]   170 64078#


exprSQ <- t(exprSQ)
exprSQ <- exprSQ[as.character(rownames(SQEIcounts)),]

geneCommon <- colnames(exprSQ)[which(is.element(as.character(colnames(exprSQ)),as.character(colnames(SQEIcounts))))]
SQEIcounts <- SQEIcounts[,geneCommon]
exprSQ <- exprSQ[,geneCommon]
rm(geneCommon)

## we check that the order of the column and rownames have the same order
stopifnot(identical(colnames(exprSQ),colnames(SQEIcounts)))
stopifnot(identical(rownames(exprSQ),rownames(SQEIcounts)))

## we calculate the intronic reads
intronicReads <- SQEIcounts - exprSQ
intronicReads[is.na(intronicReads)]=0
## dim(intronicReads)
## [1]   170 64078

summedInt <- apply(intronicReads,2,sum)
summedIntExo <- apply(SQEIcounts,2,sum)
summedExo <- apply(exprSQ,2,sum)
  
ratiosInt <- summedInt/summedIntExo
ratiosExo <- summedExo/summedIntExo

ratiosInt <- ratiosInt[-which(is.na(ratiosInt))]
ratiosExo <- ratiosExo[names(ratiosInt)]


## I think this is better to look it after RPKM convertion

jpeg("plots/intronicExonicrations.jpeg", type="cairo")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
hist(ratiosInt,col='skyblue',border=F,main= "Proportion of genic reads" ,xlab="",breaks=30)
hist(ratiosExo,add=T,col='red',border=F,xlab= "Samples",xaxt="n",breaks=30)
legend("topright", inset=c(-0.3,0), legend=c("Intronic","Exonic"), 
       fill=c("skyblue","red"), col=c("skyblue","red"), title="Reads type")
dev.off()

## Now we correct for PEER using simple quantification Exons+Introns
save(intronicReads,file="data/expr/rawCounts/genic/exprIntrons.rda")


## save the intronic expression
rm(list=ls())
setwd("/home/guelfi/eQTLPipeline")
sink("logIntronic.log")
nCores <- 7
cat(paste("Number of cores",nCores,"\n"))
library(devtools)
load_all()

## Now we correct for PEER using simple quantification Exons+Introns
load(file="data/expr/rawCounts/genic/exprIntrons.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")
intronicReads <- t(intronicReads)

## convert the genes that have NAs
intronicReads[is.na(intronicReads)]=0
table(is.na(intronicReads))
## clean negative reads
table(intronicReads[intronicReads<0])
intronicReads[intronicReads<0]=0
table(intronicReads[intronicReads<0])


## remove genes that not expressed in any gene
intronicReads <- intronicReads[rowSums(intronicReads>0)>0,]
cat("Processing PUTM \n")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]

# now we select the expression for the PUTM only samples
expr <- intronicReads[,as.character(PUTM$A.CEL_file)]
rm(intronicReads)

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content genic and gene length

#geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
#geneLength <- geneLength[as.character(rownames(expr)),c(3,1:2)]

# load the definition of the exons
exonsdef <- read.csv("data/general/exonDef.csv")

## calculation of genes only exons length


library(doParallel)
library(foreach)

detectCores()
## [1] 24

cl <- makeCluster(nCores)
clusterExport(cl,c("getIntronicRegions","defExonicRegions"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
intronicRegions <- foreach(i=1:length(rownames(expr)),.combine=rbind,.verbose=F)%dopar%getIntronicRegions(rownames(expr)[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)


write.table(data.frame(intronicRegions), file = paste0("data/general/intronicRegions.BED"), row.names = F, 
            col.names = F, quote = F)



intronDef <- read.delim("data/general/intronicRegions.BED",header=F)
colnames(intronDef) <- c("Chromosome.Name","Exon.Chr.Start..bp.","Exon.Chr.End..bp.","Ensembl.Gene.ID")

library(doParallel)
library(foreach)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(nCores)
clusterExport(cl, c("getRegionsWidth","defExonicRegions"))

registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
geneswidth <- foreach(i=1:length(rownames(expr)),.export=c("getRegionsWidth","defExonicRegions"),.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(expr)[i],intronDef,intronic=TRUE)
#geneswidth <- foreach(i=1:10,.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(expr)[i],exonsdef)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)



length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[as.character(rownames(expr))]
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

## filtering
RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std)
## write log
cat(paste("Number of Genes after filtering:",length(genesList),"\n"))
## we remove the genes that are composed only by 1 exon
genesList <- genesList[!is.na(genesList)]
expr <- expr[as.character(genesList),]

rm(RPKM.std,geneswidth)

## we have already save the intronic regions we remove the NAs
system("grep -v HG* data/general/intronicRegions.BED  | grep -v LRG* | grep -v HS* | grep -v NA | cat > data/general/intronicRegionsFiltered.BED ")

cmd <- paste0("/apps/BEDTools/2.24.0/bin/bedtools nuc -fi /home/ukbec/bowtie2Index/genome37.72.fa -bed data/general/intronicRegionsFiltered.BED > data/general/GCcontRegionsIntronic")

## calculate GC content with bedtools
system(cmd)

GCcontentTab <- read.delim("data/general/GCcontRegionsIntronic")
cat(paste("GC content saved in data/general/GCcontRegionsIntronic","\n"))

rm(cmd)

cl <- makeCluster(nCores)
clusterExport(cl, c("ratioGCcontent"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
GCcontentByGene <- foreach(i=1:length(rownames(expr)),.combine=rbind,.verbose=F)%dopar%ratioGCcontent(rownames(expr)[i],GCcontentTab)
end <- Sys.time()
stopCluster(cl)
cat(paste("GC content ratios calculated in ",end-start,"\n"))
rm(cl,end,start,GCcontentTab,exonsdef)

GCcontent <- GCcontentByGene
rownames(GCcontent) <- GCcontent[,1]
GCcontent <- as.data.frame(GCcontent[,-1]) 
length <- length[as.character(rownames(GCcontent))]  
stopifnot(identical(names(length),rownames(GCcontent)))
GCcontent <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
colnames(GCcontent) <- c("GCcontent","length")
# we update gene expression with the filtered genes
rm(length,geneLength,genesList)
GCcontent$GCcontent <- as.numeric(GCcontent$GCcontent)

library(cqn)
library(scales)

## CQN
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),rownames(GCcontent)))              

cat("Conditional quantile normalisation \n")
my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent,sizeFactors=librarySize, verbose = TRUE)


png(paste0("plots/introns/CQNPUTM.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset


cat(paste("Number of Genes and samples",dim(RPKM.cqn),"\n"))
# 25985 genic regions kept
# length(grep("DER",rownames(RPKM.cqn)))
# 42577 intergenic regiond kept

# save results
cat("Saving the the RPKM CQN normalised in data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
##save(RPKM.cqn,file="data/expr/normalisedCounts/SQPRKM.cqn.rda",compress="bzip2")

PUTM$U.Region_simplified <- NULL
covs <- PUTM 
rownames(covs) <- covs$A.CEL_file
#convert the female and male info in numeric
covs[covs=="M"]=0
covs[covs=="F"]=1
covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
covs[,5] <- as.numeric(covs[,5])
covs[,6] <- as.numeric(covs[,6])
lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
rownames(lanes) <- gsub("CEL","",rownames(lanes))
covs <- cbind(covs,librarySize[as.character(rownames(covs))])
covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                    "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")

save(RPKM.cqn,PUTM,covs,file="data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")



### Residual correction ###

rm(list=ls())

load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")

##doSwamp(RPKM.cqn,covs)

PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                               "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda")

##doSwamp(resids,covs)



############
### SNIG ###
############



## save the intronic expression
rm(list=ls())
setwd("/home/guelfi/eQTLPipeline")
sink("logIntronic.log")
nCores <- 7
cat(paste("Number of cores",nCores,"\n"))
library(devtools)
load_all()

## Now we correct for PEER using simple quantification Exons+Introns
load(file="data/expr/rawCounts/genic/exprIntrons.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")
intronicReads <- t(intronicReads)

## convert the genes that have NAs
intronicReads[is.na(intronicReads)]=0
table(is.na(intronicReads))
## clean negative reads
table(intronicReads[intronicReads<0])
intronicReads[intronicReads<0]=0
table(intronicReads[intronicReads<0])


## remove genes that not expressed in any gene
intronicReads <- intronicReads[rowSums(intronicReads>0)>0,]
cat("Processing SNIG \n")
SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),]

# now we select the expression for the SNIG only samples
expr <- intronicReads[,as.character(SNIG$A.CEL_file)]
rm(intronicReads)

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(SNIG$A.CEL_file),]
names(librarySize) <- as.character(SNIG$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content genic and gene length

#geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
#geneLength <- geneLength[as.character(rownames(expr)),c(3,1:2)]

# load the definition of the exons
exonsdef <- read.csv("data/general/exonDef.csv")

## calculation of genes only exons length


library(doParallel)
library(foreach)

detectCores()
## [1] 24

cl <- makeCluster(nCores)
clusterExport(cl,c("getIntronicRegions","defExonicRegions"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
intronicRegions <- foreach(i=1:length(rownames(expr)),.combine=rbind,.verbose=F)%dopar%getIntronicRegions(rownames(expr)[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)


write.table(data.frame(intronicRegions), file = paste0("data/general/intronicRegions.BED"), row.names = F, 
            col.names = F, quote = F)




intronDef <- read.delim("data/general/intronicRegions.BED",header=F)
colnames(intronDef) <- c("Chromosome.Name","Exon.Chr.Start..bp.","Exon.Chr.End..bp.","Ensembl.Gene.ID")

library(doParallel)
library(foreach)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(nCores)
clusterExport(cl, c("getRegionsWidth","defExonicRegions"))

registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
geneswidth <- foreach(i=1:length(rownames(expr)),.export=c("getRegionsWidth","defExonicRegions"),.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(expr)[i],intronDef,intronic=TRUE)
#geneswidth <- foreach(i=1:10,.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(expr)[i],exonsdef)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)



length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[as.character(rownames(expr))]
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

## filtering
RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std)
## write log
cat(paste("Number of Genes after filtering:",length(genesList),"\n"))
## we remove the genes that are composed only by 1 exon
genesList <- genesList[!is.na(genesList)]
expr <- expr[as.character(genesList),]

rm(RPKM.std,geneswidth)

## we have already save the intronic regions we remove the NAs
system("grep -v HG* data/general/intronicRegions.BED  | grep -v LRG* | grep -v HS* | grep -v NA | cat > data/general/intronicRegionsFiltered.BED ")

cmd <- paste0("/apps/BEDTools/2.24.0/bin/bedtools nuc -fi /home/ukbec/bowtie2Index/genome37.72.fa -bed data/general/intronicRegionsFiltered.BED > data/general/GCcontRegionsIntronic")

## calculate GC content with bedtools
system(cmd)

GCcontentTab <- read.delim("data/general/GCcontRegionsIntronic")
cat(paste("GC content saved in data/general/GCcontRegionsIntronic","\n"))

rm(cmd)

cl <- makeCluster(nCores)
clusterExport(cl, c("ratioGCcontent"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
GCcontentByGene <- foreach(i=1:length(rownames(expr)),.combine=rbind,.verbose=F)%dopar%ratioGCcontent(rownames(expr)[i],GCcontentTab)
end <- Sys.time()
stopCluster(cl)
cat(paste("GC content ratios calculated in ",end-start,"\n"))
rm(cl,end,start,GCcontentTab,exonsdef)

GCcontent <- GCcontentByGene
rownames(GCcontent) <- GCcontent[,1]
GCcontent <- as.data.frame(GCcontent[,-1]) 
length <- length[as.character(rownames(GCcontent))]  
stopifnot(identical(names(length),rownames(GCcontent)))
GCcontent <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
colnames(GCcontent) <- c("GCcontent","length")
# we update gene expression with the filtered genes
rm(length,geneLength,genesList)
GCcontent$GCcontent <- as.numeric(GCcontent$GCcontent)

library(cqn)
library(scales)

## CQN
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),rownames(GCcontent)))              

cat("Conditional quantile normalisation \n")
my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent,sizeFactors=librarySize, verbose = TRUE)


png(paste0("plots/introns/CQNSNIG.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset


cat(paste("Number of Genes and samples",dim(RPKM.cqn),"\n"))
# 25985 genic regions kept
# length(grep("DER",rownames(RPKM.cqn)))
# 42577 intergenic regiond kept

# save results
cat("Saving the the RPKM CQN normalised in data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.SNIG")
##save(RPKM.cqn,file="data/expr/normalisedCounts/SQPRKM.cqn.rda",compress="bzip2")

SNIG$U.Region_simplified <- NULL
covs <- SNIG 
rownames(covs) <- covs$A.CEL_file
#convert the female and male info in numeric
covs[covs=="M"]=0
covs[covs=="F"]=1
covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
covs[,5] <- as.numeric(covs[,5])
covs[,6] <- as.numeric(covs[,6])
lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
rownames(lanes) <- gsub("CEL","",rownames(lanes))
covs <- cbind(covs,librarySize[as.character(rownames(covs))])
covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                    "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")

save(RPKM.cqn,SNIG,covs,file="data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.SNIG")


### PCA after normalisation ###

load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.SNIG")
RPKM.cqn.SNIG <- RPKM.cqn
load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
RPKM.cqn.PUTM <- RPKM.cqn

rm(RPKM.cqn)

comJunc <- intersect(rownames(RPKM.cqn.SNIG),rownames(RPKM.cqn.PUTM))
length(comJunc)
RPKM.cqn <- cbind(RPKM.cqn.SNIG[as.character(comJunc),],RPKM.cqn.PUTM[as.character(comJunc),])

PCAres<- prcomp(t(RPKM.cqn))
par(mfrow=c(1,1))

plot(PCAres, main="PCA axes gene-intronic(PUTM + SNIG)")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue (gene-intronic)",xlab="PC1",ylab="PC2",ylim=c(-350,250),xlim=c(-500,850))

points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    


rm(PUTM,SNIG,PCAres)


### Residual correction ###

rm(list=ls())

load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.SNIG")

##doSwamp(RPKM.cqn,covs)

PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                               "data/expr/normalisedCounts/genic/geneIntronic/resids.SNIG.rda")



setwd("/home/guelfi/eQTLPipeline")


### Run the eQTL analysis
writeSH(nameSH="runCisEQTL.sh",logName="runCisEQTL",
        cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/parRunCiseQTLIntronic.R"),numThreads=8)


system("qsub runCisEQTL.sh")


### Run the Sentinalisation
writeSH(nameSH="LDsentinalisation.sh",logName="LDsentinalisation",
        cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/sentiIntronic.R"),numThreads=8)

system("qsub LDsentinalisation.sh")





eQTLPUTM <- read.delim("data/results/finaleQTLs/intronic.PUTM.txt",sep=" ")
table(eQTLPUTM$myFDR<0.05)
#     FALSE  TRUE
#     494  1795
table(eQTLPUTM$myFDR<0.01)
#     FALSE  TRUE
#      1222  1067


eQTLSNIG <- read.delim("data/results/finaleQTLs/intronic.SNIG.txt",sep=" ")
table(eQTLSNIG$myFDR<0.05)
#     FALSE  TRUE
#     430  1015
table(eQTLSNIG$myFDR<0.01)
#     FALSE  TRUE
#     886   559





###########################
####    Annotation ########
###########################

### PUTM ###


library("biomaRt")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","gene_biotype","description"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=eQTLPUTM$gene, mart=ensembl)



rownames(geneNames)<- geneNames$ensembl_gene_id
geneNames$ensembl_gene_id <- NULL


eQTLPUTM <- cbind(eQTLPUTM,geneNames[as.character(eQTLPUTM$gene),])

save(eQTLPUTM,file="data/results/finaleQTLs/intronic.Ann.PUTM.rda")
write.csv(eQTLPUTM,file="data/results/finaleQTLs/intronic.Ann.PUTM.csv",row.names=F)    

### SNIG ###

rm(geneNames)

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","gene_biotype","description"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=eQTLSNIG$gene, mart=ensembl)



rownames(geneNames)<- geneNames$ensembl_gene_id
geneNames$ensembl_gene_id <- NULL


eQTLSNIG <- cbind(eQTLSNIG,geneNames[as.character(eQTLSNIG$gene),])

save(eQTLSNIG,file="data/results/finaleQTLs/intronic.Ann.SNIG.rda")
write.csv(eQTLSNIG,file="data/results/finaleQTLs/intronic.Ann.SNIG.csv",row.names=F)    


#########################
### Annotation by SNP ###
#########################


rm(list=ls())

load(file="data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTM <-  annSNP(eQTLPUTM)

save(eQTLPUTM,file="data/results/finaleQTLs/intronic.Ann.PUTM.rda")
write.csv(eQTLPUTM,file="data/results/finaleQTLs/intronic.Ann.PUTM.csv",row.names=F)    

load(file="data/results/finaleQTLs/intronic.Ann.SNIG.rda")
eQTLSNIG <-  annSNP(eQTLSNIG)

save(eQTLSNIG,file="data/results/finaleQTLs/intronic.Ann.SNIG.rda")
write.csv(eQTLSNIG,file="data/results/finaleQTLs/intronic.Ann.SNIG.csv",row.names=F)    


########################
### Annotate the TSS ###
########################

library("biomaRt")


load(file="data/results/finaleQTLs/intronic.Ann.PUTM.rda")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=eQTLPUTM$gene, mart=ensembl)


TSS <- sapply(eQTLPUTM$gene, function(x){getTSS(x,geneNames)})

eQTLPUTM <- cbind(eQTLPUTM,TSS)

save(eQTLPUTM,file="data/results/finaleQTLs/intronic.Ann.PUTM.rda")
write.csv(eQTLPUTM,file="data/results/finaleQTLs/intronic.Ann.PUTM.rda",row.names=F)    
rm(eQTLPUTM,TSS,geneNames)

### SNIG

load(file="data/results/finaleQTLs/intronic.Ann.SNIG.rda")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=eQTLSNIG$gene, mart=ensembl)

TSS <- sapply(eQTLSNIG$gene, function(x){getTSS(x,geneNames)})

eQTLSNIG <- cbind(eQTLSNIG,TSS)

save(eQTLSNIG,file="data/results/finaleQTLs/intronic.Ann.SNIG.rda")
write.csv(eQTLSNIG,file="data/results/finaleQTLs/intronic.Ann.SNIG.rda",row.names=F)    









 
