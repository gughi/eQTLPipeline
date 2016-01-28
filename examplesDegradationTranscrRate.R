## example for Mina

par(mfrow=c(1,2))

gene <- "ENSG00000267053"
snp <- "chr19:36875996:TG_T"
##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
#exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="PUTM", main=paste("exonic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )



gene <- "ENSG00000267053"
snp <- "chr19:36875996:TG_T"
##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

#exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="PUTM", main=paste("intronic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )




#####################################
### chr10:51801386 ENSG00000099290 ##
#####################################

## example for Mina

par(mfrow=c(1,2))

gene <- "ENSG00000099290"
snp <- "chr10:51801386"
##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
#exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="PUTM", main=paste("exonic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )



##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

#exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="PUTM", main=paste("intronic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )







## example for Mina

par(mfrow=c(1,2))

gene <- "ENSG00000100376"
snp <- "chr22:45731730"
##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
#exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="PUTM", main=paste("exonic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )



##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

#exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="PUTM", main=paste("intronic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )










## example for JUAN SP1
# library(devtools)
# load_all()
par(mfrow=c(1,2))

gene <- "ENSG00000185591"
snp <- "chr12:53513165"
##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.SNIG.rda"
##exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="SNIG", main=paste("exonic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )



##snp <- "chr6:29955809"
# dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
# infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"

dosageFile <- "/home/adai/genotyped/imputed_v3/imputed.dosage"
infoFile <- "/home/adai/genotyped/imputed_v3/imputed.info"

#exprFile <- "data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda"
exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.SNIG.rda"

load("data/general/sampleInfo.rda")
eigenFile <- "/home/seb/plinkOutput/eigenvec"
title <- paste0(gene ,"-",snp)

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

plot(round(PUTMsnp),PUTMtmp,xaxt="n",ylab="expression",xlab="SNIG", main=paste("intronic\n",title),cex.main=0.8)
abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
fit <- coef( summary(glm( PUTMtmp ~ PUTMsnp + my.cov0[1,]+my.cov0[2,]+my.cov0[3,]) ))

mtext(paste("pval:",fit["PUTMsnp", "Pr(>|t|)"],"beta:",fit["PUTMsnp", "Estimate"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals

mtext( side=1, paste0(info["Al2"], info["Al2"]), at=(0), cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al2"], info["Al1"]), at=(1),   cex=0.6, line=-0.25 )
mtext( side=1, paste0(info["Al1"], info["Al1"]), at=(2), cex=0.6, line=-0.25 )


