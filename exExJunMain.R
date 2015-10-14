
# Apollo
##setwd("/home/guelfi/eQTLPipeline")
# Caprica
setwd("/home/seb/projectsR/eQTLPipeline")
sink("logExonExonJunctions.log")
nCores <- 15
cat(paste("Number of cores",nCores,"\n"))
library(devtools)
load_all()

# Now we correct for PEER using simple quantification Exons+Introns



# expr <- read.delim("data/expr/rawCounts/genic/exonExonJunctions")
# ## first four columns have the exon1ID - exon2ID - chr and TSS
# 
# ## Loading information of exon -exon junctions
# # chr - start - end - ID
# map <- read.delim(pipe("grep GB data/general/exongrps.log | cut -f2-5 -d' '"),sep=" ",header=F)
# # Number of exons - length - max exon length If it's a complicated groupr and whether there are overlapping exons
# mapTmp <- read.delim(pipe("grep GI data/general/exongrps.log | cut -f2-6 -d' '"),sep=" ",header=F)
# map <- cbind(map,mapTmp)
# 
# ## load the map of the ID with the exon IDs
# mapTmp <- read.delim(pipe("grep GE data/general/exongrps.log"),sep=" ",header=F,fill=T , col.names=paste("col", 1:64, sep="_"))
# ##mapTmp[1,!is.na(mapTmp[1,])]
# 
# y <-1:nrow(mapTmp)
# IDs <- unlist(sapply(y, function(x){ rep(x,length(mapTmp[x,2:64][!is.na(mapTmp[x,2:64])]))}))
# 
# mapExon <- read.delim(pipe("grep -w E data/general/exongrps.log"),sep=" ",header=F)
# mapExon <- cbind(mapExon,IDs)
# ## columns check the Altrans manual if something unclear
# mapExon$V1 <- NULL
# 
# 
# ## clean the data
# 
# colnames(mapExon) <- c("chr","start", "end", "exonID","geneID", "transID", "unifiedExon",
#                      "differentStrand","strand","length","UnRegStartEnd",
#                      "groupID")
# 
# colnames(map) <- c("chr","start", "end", "groupID","noExons", "length",
#                    "maxExonLength","probleGroup","nonOverlaExon")
# 
# map$noExons <- gsub("NoExons:","",map$noExons)
# map$length <- gsub("Length:","",map$length)
# map$maxExonLength <- gsub("MaxExonLength:","",map$maxExonLength)
# save(map,mapExon,expr,file="data/expr/rawCounts/fullExExJun.rda")

load("data/expr/rawCounts/genic/fullExExJun.rda")

# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")
## reformat the sample names to match exon-exon Junction expression names

# investigate the PC1
# PCA PC1 vs PC2 
PCAres<- prcomp(t(expr[,5:ncol(expr)]))
par(mfrow=c(1,1))
PUTM <- sampleInfo[sampleInfo[, 6] == "PUTM",]
SNIG <- sampleInfo[sampleInfo[, 6] == "SNIG",]

plot(PCAres, main="PCA axis exon-exon juctions (PUTM + SNIG)")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue",xlab="PC1",ylab="PC2" )

exprSamNam <- sapply(PUTM$A.CEL_file,function(x){
  paste(unlist(strsplit(x,"_")),collapse="_")})

points(PCAres$x[paste0("Sample_",exprSamNam),1],
       PCAres$x[paste0("Sample_",exprSamNam),2],col="red")

exprSamNam <- sapply(SNIG$A.CEL_file,function(x){
  paste(unlist(strsplit(x,"_")),collapse="_")})

points(PCAres$x[paste0("Sample_",exprSamNam),1],
       PCAres$x[paste0("Sample_",exprSamNam),2],col="blue")

legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    

rm(PUTM,SNIG,PCAres)

doSwamp(t(expr[,5:ncol(expr)]),covs=sampleInfo)
## we calculate the GC content for all the junctions

write.table(data.frame(mapExon[,1:4]), file = paste0("data/general/regionJunction.BED"), row.names = F, 
            col.names = F, quote = F,sep="\t")

## we filter things that not match with the fasta file

system("grep -v HG* data/general/regionJunction.BED  | grep -v LRG* | grep -v HS* | cat > data/general/regionJunctionFiltered.BED")

cmd <- paste0("/apps/BEDTools/2.24.0/bin/bedtools nuc -fi /home/ukbec/bowtie2Index/genome37.72.fa -bed data/general/regionJunctionFiltered.BED > data/general/GCcontRegionJunctions; rm data/general/regionJunctionFiltered.BED ;rm data/general/regionJunction.BED")

## calculate GC content with bedtools
system(cmd)
rm(cmd,exprSamNam)

cat(paste("GC content saved in data/general/GCcontRegionJunctions","\n"))
GCcontentTab <- read.delim("data/general/GCcontRegionJunctions")



# now we select the expression for the PUTM only samples
exprAll <- expr
colnames(expr) <- paste0(gsub("Sample_","",colnames(expr)),"_")

# now we select the expression for the PUTM only samples
expr <- expr[,as.character(PUTM$A.CEL_file)]

## detectCores()
## [1] 24
librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content genic and gene length

juncdef <- exprAll[,1:4]

## calculation of genes only exons length
# THIS WAS CALCULATED IN CAPRICA BECAUSE didn't manage to do it in Apollo
#
# library(doParallel)
# library(foreach)
# 
# detectCores()
# ## [1] 24
# # create the cluster with the functions needed to run
# cl <- makeCluster(20)
# clusterExport(cl, c("lengthJunction"))
# 
# registerDoParallel(cl)
# getDoParWorkers()
# 
# start <- Sys.time()
# system.time(length <- foreach(i=1:nrow(juncdef[,1:2]),.combine=c,.verbose=F)%dopar%lengthJunction(juncdef[i,1:2],mapExon))
# #geneswidth <- foreach(i=1:10,.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(expr)[i],exonsdef)
# ##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
# end <- Sys.time()
# end-start
# stopCluster(cl)
# rm(cl,end,start)
# 
# 
# length <- apply(juncdef[,1:2],1,function(x) lengthJunction(x,mapExon))


load("data/general/lengthExExJun.rda")
juncdef <- cbind(juncdef,length)
IDs <- do.call(paste, c(juncdef[,1:2],sep="_"))
rownames(juncdef) <- IDs
rownames(expr) <- IDs
names(length) <- IDs

stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              


RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

## BEFORE FILTERING WE HAD 551102 exon-exon junctions

## filtering
RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std)

# ALTER FILTERING 96738 RETAINED exon-exon junctions
## write log
cat(paste("Number of exon-exon after filtering:",length(genesList),"\n"))
rm(RPKM.std)
expr <- expr[as.character(genesList),]

library(doParallel)
library(foreach)

cl <- makeCluster(15)
clusterExport(cl, c("ratioGCcontentExonExon"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
GCcontentByExEx <- foreach(i=1:length(rownames(expr)),.combine=rbind,.verbose=F)%dopar%ratioGCcontentExonExon(rownames(expr)[i],GCcontentTab)
end <- Sys.time()
stopCluster(cl)
cat(paste("GC content ratios calculated in ",end-start,"\n"))


GCcontent <- GCcontentByExEx
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

png(paste0("plots/exonExonJunc/CQNPUTM.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset

cat(paste("Number of Genes and samples",dim(RPKM.cqn),"\n"))

# save results
cat("Saving the the RPKM CQN normalised in data/expr/normalisedCounts/RPKM.cqn.PUTM \n")
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

save(RPKM.cqn,PUTM,covs,file="data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.PUTM")


### Residual correction ###

rm(list=ls())

load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.PUTM")

##doSwamp(RPKM.cqn,covs)

PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                               "data/expr/normalisedCounts/genic/exonExonJunc/resids.PUTM.rda")





##############
### SNIG #####
##############



# now we select the expression for the SNIG only samples
expr <- exprAll
colnames(expr) <- paste0(gsub("Sample_","",colnames(expr)),"_")

# now we select the expression for the SNIG only samples
expr <- expr[,as.character(SNIG$A.CEL_file)]

## detectCores()
## [1] 24
librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(SNIG$A.CEL_file),]
names(librarySize) <- as.character(SNIG$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content genic and gene length

juncdef <- exprAll[,1:4]

load("data/general/lengthExExJun.rda")
juncdef <- cbind(juncdef,length)
IDs <- do.call(paste, c(juncdef[,1:2],sep="_"))
rownames(juncdef) <- IDs
rownames(expr) <- IDs
names(length) <- IDs

stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              


RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

## BEFORE FILTERING WE HAD 551102 exon-exon junctions

## filtering
RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std)

# ALTER FILTERING 88936 RETAINED exon-exon junctions
## write log
cat(paste("Number of exon-exon after filtering:",length(genesList),"\n"))
rm(RPKM.std)
expr <- expr[as.character(genesList),]

library(doParallel)
library(foreach)

cl <- makeCluster(15)
clusterExport(cl, c("ratioGCcontentExonExon"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
GCcontentByExEx <- foreach(i=1:length(rownames(expr)),.combine=rbind,.verbose=F)%dopar%ratioGCcontentExonExon(rownames(expr)[i],GCcontentTab)
end <- Sys.time()
stopCluster(cl)
cat(paste("GC content ratios calculated in ",end-start,"\n"))


GCcontent <- GCcontentByExEx
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

png(paste0("plots/exonExonJunc/CQNSNIG.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset

cat(paste("Number of Genes and samples",dim(RPKM.cqn),"\n"))

# save results
cat("Saving the the RPKM CQN normalised in data/expr/normalisedCounts/RPKM.cqn.SNIG \n")
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

save(RPKM.cqn,SNIG,covs,file="data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.SNIG")


### Residual correction ###

rm(list=ls())

load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.SNIG")

##doSwamp(RPKM.cqn,covs)

PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                               "data/expr/normalisedCounts/genic/exonExonJunc/resids.SNIG.rda")




## check PCA before correcting

load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.SNIG")
RPKM.cqn.SNIG <- RPKM.cqn
load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.PUTM")
RPKM.cqn.PUTM <- RPKM.cqn

rm(RPKM.cqn)

comJunc <- intersect(rownames(RPKM.cqn.SNIG),rownames(RPKM.cqn.PUTM))
length(comJunc)
RPKM.cqn <- cbind(RPKM.cqn.SNIG[as.character(comJunc),],RPKM.cqn.PUTM[as.character(comJunc),])

PCAres<- prcomp(t(RPKM.cqn))
par(mfrow=c(1,1))
PUTM <- sampleInfo[sampleInfo[, 6] == "PUTM",]
SNIG <- sampleInfo[sampleInfo[, 6] == "SNIG",]

plot(PCAres, main="PCA axis exon-exon juctions (PUTM + SNIG)")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue",xlab="PC1",ylab="PC2" )

points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    

rm(PUTM,SNIG,PCAres)





 
  
  
  
  