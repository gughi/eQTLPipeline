## Now we check the intergenic regions

## So we subdivide the intergenic regions with eQTls and the intergenic 

## subdivide the groups
load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

## 5%FDR
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.05),1:7]
head(eQTLPUTM)

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
RPKM.cqn.PUTM <- RPKM.cqn
defReg.PUTM <- starStopReg
rm(PUTM,RPKM.cqn,covs)

eQTLPUTM <- cbind(eQTLPUTM,starStopReg[as.character(eQTLPUTM[,2]),])

load("data/results/finaleQTLs/intergenic.Ann.SNIG.rda")

## 5%FDR
eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.05),1:7]
head(eQTLSNIG)

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.SNIG")
RPKM.cqn.SNIG <- RPKM.cqn
defReg.SNIG <- starStopReg
rm(SNIG,RPKM.cqn,covs)

eQTLSNIG <- cbind(eQTLSNIG,starStopReg[as.character(eQTLSNIG[,2]),])
rm(starStopReg)


tmpRegiosPUTM <- GRanges(defReg.PUTM$chr, IRanges(defReg.PUTM$start, defReg.PUTM$end),)
tmpRegiosSNIG <- GRanges(defReg.SNIG$chr, IRanges(defReg.SNIG$start, defReg.SNIG$end))
names(tmpRegiosPUTM) <- rownames(defReg.PUTM)
names(tmpRegiosSNIG) <- rownames(defReg.SNIG)
head(tmpRegiosPUTM)
head(tmpRegiosSNIG)
## select the transcribed regions identified in the data

overlap <- findOverlaps(tmpRegiosSNIG,tmpRegiosPUTM)

overlap <- as.data.frame(overlap)
## we create the mapping file

idxSNIG <- names(tmpRegiosSNIG[overlap[,1]])
idxPUTM <- names(tmpRegiosPUTM[overlap[,2]])

overlap <- cbind(overlap,names(tmpRegiosSNIG[overlap[,1]]),names(tmpRegiosPUTM[overlap[,2]]))


mergedIntergenic <- cbind(RPKM.cqn.PUTM[idxPUTM,],RPKM.cqn.SNIG[idxSNIG,])


idx.eQTL.PUTM <- intersect(idxPUTM,as.character(eQTLPUTM$gene))
length(idx.eQTL.PUTM)                           
## only 1222 eQTL regions were identify as transcribed regions in the SNIG dataset


idx.eQTL.SNIG <- intersect(idxSNIG,as.character(eQTLSNIG$gene))
length(idx.eQTL.SNIG)                           
## only 681 eQTL regions were identify as transcribed regions in the SNIG dataset


idxEQTL <- overlap[which(as.character(overlap[,4]) %in% as.character(idx.eQTL.PUTM)),]

idxEQTL <- rbind(idxEQTL,overlap[which(as.character(overlap[,3]) %in% as.character(idx.eQTL.SNIG)),])

## remove the eQTLs that overlap
idxEQTL <- unique(idxEQTL)


exprEQTL.int <- mergedIntergenic[as.character(idxEQTL[,4]),]
exprNonEQTL.int <- mergedIntergenic[-which(rownames(mergedIntergenic) %in% as.character(idxEQTL[,4])),]

hist(exprEQTL.int,col='skyblue',border=F,main= "expression eQTL vs non-eQTL",freq=F
     ,breaks = 40,xlab="RPKM.cqn")
hist(exprNonEQTL.int,add=T,col=scales::alpha('red',.5),border=F,freq=F,breaks=40)
lines(density(exprEQTL.int, adjust = 2), col = "skyblue",lwd=3)
lines(density(exprNonEQTL.int, adjust = 2), col = "red",lwd=3)
legend("topright",c("eQTL","non-eQTL"),col=c('skyblue','red'),pch=15)





load("data/general/sampleInfo.rda")

PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),]


## histogram of the expression divided by tissue and eQTL and non-eQTL regions

##PUTM
exprEQTL.int <- mergedIntergenic[as.character(idxEQTL[,4]),as.character(PUTM$A.CEL_file)]
exprNonEQTL.int <- mergedIntergenic[-which(rownames(mergedIntergenic) %in% as.character(idxEQTL[,4])),as.character(PUTM$A.CEL_file)]

hist(exprEQTL.int,col='skyblue',border=F,main= "expression eQTL vs non-eQTL (PUTM)",freq=F
     ,breaks = 40,xlab="RPKM.cqn")
hist(exprNonEQTL.int,add=T,col=scales::alpha('red',.5),border=F,freq=F,breaks=40)
lines(density(exprEQTL.int, adjust = 2), col = "skyblue",lwd=3)
lines(density(exprNonEQTL.int, adjust = 2), col = "red",lwd=3)
legend("topright",c("eQTL","non-eQTL"),col=c('skyblue','red'),pch=15)

## SNIG

exprEQTL.int <- mergedIntergenic[as.character(idxEQTL[,4]),as.character(SNIG$A.CEL_file)]
exprNonEQTL.int <- mergedIntergenic[-which(rownames(mergedIntergenic) %in% as.character(idxEQTL[,4])),as.character(SNIG$A.CEL_file)]

hist(exprEQTL.int,col='skyblue',border=F,main= "expression eQTL vs non-eQTL (SNIG)",freq=F
     ,breaks = 40,xlab="RPKM.cqn")
hist(exprNonEQTL.int,add=T,col=scales::alpha('red',.5),border=F,freq=F,breaks=40)
lines(density(exprEQTL.int, adjust = 2), col = "skyblue",lwd=3)
lines(density(exprNonEQTL.int, adjust = 2), col = "red",lwd=3)
legend("topright",c("eQTL","non-eQTL"),col=c('skyblue','red'),pch=15)






PCAres<- prcomp(t(exprNonEQTL.int),scale=TRUE)
par(mfrow=c(1,1))

plot(PCAres, main="PCA axes intergenic regions")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue (intergenic quantification) non eQTL regions",xlab="PC1",ylab="PC2",ylim=c(-350,250),xlim=c(-500,850))

points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    


PCAres<- prcomp(t(exprEQTL.int),scale=TRUE)
par(mfrow=c(1,1))

plot(PCAres, main="PCA axes intergenic regions")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue (intergenic quantification) eQTL regions",xlab="PC1",ylab="PC2",ylim=c(-350,250),xlim=c(-500,850))

points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    
rm(exprEQTL.int,exprNonEQTL.int,mergedIntergenic,RPKM.cqn.PUTM,RPKM.cqn.SNIG,PCAres)


load("data/expr/rawCounts/intergenic/transcribedRegions/SNIG/allChromosomes.rda")

exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
colnames(exprIntergenic) <- c("chr","start","end","width",colnames(coverage))
rm(annotation,coverage)    


## select regions with length => 100bp
exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
# We don't do filtering since it doesn't make any since for DERFINDER; 
# regions are detected based on the expression
# creation of identifiers for the regions
IDs <- paste0("DER",c(1:nrow(exprIntergenic)))
rownames(exprIntergenic) <- IDs
rm(IDs)



rawExpr.SNIG <- exprIntergenic[as.character(overlap[,3]),]
rm(exprIntergenic)

load("data/expr/rawCounts/intergenic/transcribedRegions/PUTM/allChromosomes.rda")


exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
colnames(exprIntergenic) <- c("chr","start","end","width",colnames(coverage))
rm(annotation,coverage)    


## select regions with length => 100bp
exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
# We don't do filtering since it doesn't make any since for DERFINDER; 
# regions are detected based on the expression
# creation of identifiers for the regions
IDs <- paste0("DER",c(1:nrow(exprIntergenic)))
rownames(exprIntergenic) <- IDs
rm(IDs)

rawExpr.PUTM <- exprIntergenic[as.character(overlap[,4]),]
rm(exprIntergenic)

rawIntergenic <- cbind(rawExpr.PUTM[,-c(1:4)],rawExpr.SNIG[,-c(1:4)])


exprEQTL.int <- rawIntergenic[as.character(idxEQTL[,4]),as.character(c(PUTM$A.CEL_file,SNIG$A.CEL_file))]
exprNonEQTL.int <- rawIntergenic[-which(rownames(rawIntergenic) %in% as.character(idxEQTL[,4])),as.character(c(PUTM$A.CEL_file,SNIG$A.CEL_file))]


PCAres<- prcomp(t(exprEQTL.int),scale=TRUE)
par(mfrow=c(1,1))

plot(PCAres, main="PCA axes intergenic regions")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue (intergenic quantification) eQTL regions",xlab="PC1",ylab="PC2",ylim=c(-300,100),xlim=c(-200,400))

points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    

PCAres<- prcomp(t(exprNonEQTL.int),scale=TRUE)
par(mfrow=c(1,1))

plot(PCAres, main="PCA axes intergenic regions")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue (intergenic quantification) non eQTL regions",xlab="PC1",ylab="PC2",ylim=c(-300,100),xlim=c(-200,400))

points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    
rm(exprEQTL.int,exprNonEQTL.int,mergedIntergenic,RPKM.cqn.PUTM,RPKM.cqn.SNIG,PCAres)

