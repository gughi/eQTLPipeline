
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")

eQTLPUTM <- eQTLPUTM[which(as.numeric(eQTLPUTM$myFDR) < 0.01),]
eQTLSNIG <- eQTLSNIG[which(as.numeric(eQTLSNIG$myFDR) < 0.01),]

library(devtools)
load_all()

infotmp <- read.table.rows(filename="/home/adai/genotyped/imputed_v3/imputed.info",
                           gsub("chr","",eQTLPUTM$snps),sep=" ")

variantsPUTM <- infotmp[as.character(gsub("chr","",eQTLPUTM$snps)),]


tmpVariansPUTM <- as.data.frame(cbind(as.character(eQTLPUTM$rs),as.character(gsub("chr","",eQTLPUTM$snps))))
head(tmpVariansPUTM)

tmpVariansPUTM$chr <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,2]),":")
                                    ,function(x){c(x[1])}))

tmpVariansPUTM$BP <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,2]),":")
                                   ,function(x){c(x[2])}))

tmpVariansPUTM$tmp <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,2]),":")
                                    ,function(x){c(x[3])}))

tmpVariansPUTM$Allele1 <- as.character(variantsPUTM$Al1)

tmpVariansPUTM$Allele2 <- as.character(variantsPUTM$Al2)

for (i in 1:nrow(tmpVariansPUTM))
{
  if(!is.na(tmpVariansPUTM$tmp[i]))
  {
    tmp <- unlist(strsplit(as.character(tmpVariansPUTM$tmp[i]),"_"))
    
    tmpVariansPUTM$Allele1[i] <- tmp[1]
    tmpVariansPUTM$Allele2[i] <- tmp[2]
  }
}


tmpVariansPUTM$Name <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,1]),";")
                                     ,function(x){c(x[1])}))
head(tmpVariansPUTM)

head(variantsPUTM)

tmpVariansPUTM <- tmpVariansPUTM[,c(8,3,4,6,7)]

colnames(tmpVariansPUTM) <- c("Name","CHR","BP","Allele1","Allele2")

write.csv(tmpVariansPUTM,file="variants.PUTM.csv",row.names=F)

head(tmpVariansPUTM)
### SNIG

infotmp <- read.table.rows(filename="/home/adai/genotyped/imputed_v3/imputed.info",
                           gsub("chr","",eQTLSNIG$snps),sep=" ")

variantsSNIG <- infotmp[as.character(gsub("chr","",eQTLSNIG$snps)),]


tmpVariansSNIG <- as.data.frame(cbind(as.character(eQTLSNIG$rs),as.character(gsub("chr","",eQTLSNIG$snps))))
head(tmpVariansSNIG)

tmpVariansSNIG$chr <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,2]),":")
                                    ,function(x){c(x[1])}))

tmpVariansSNIG$BP <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,2]),":")
                                   ,function(x){c(x[2])}))

tmpVariansSNIG$tmp <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,2]),":")
                                    ,function(x){c(x[3])}))

tmpVariansSNIG$Allele1 <- as.character(variantsSNIG$Al1)

tmpVariansSNIG$Allele2 <- as.character(variantsSNIG$Al2)

for (i in 1:nrow(tmpVariansSNIG))
{
  if(!is.na(tmpVariansSNIG$tmp[i]))
  {
    tmp <- unlist(strsplit(as.character(tmpVariansSNIG$tmp[i]),"_"))
    
    tmpVariansSNIG$Allele1[i] <- tmp[1]
    tmpVariansSNIG$Allele2[i] <- tmp[2]
  }
}


tmpVariansSNIG$Name <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,1]),";")
                                     ,function(x){c(x[1])}))
head(tmpVariansSNIG)

head(variantsSNIG)

tmpVariansSNIG <- tmpVariansSNIG[,c(8,3,4,6,7)]

colnames(tmpVariansSNIG) <- c("Name","CHR","BP","Allele1","Allele2")

write.csv(tmpVariansSNIG,file="variants.SNIG.csv",row.names=F)

rm(listìls())

load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")

eQTLPUTM <- eQTLPUTM[which(as.numeric(eQTLPUTM$myFDR) < 0.01),]
eQTLSNIG <- eQTLSNIG[which(as.numeric(eQTLSNIG$myFDR) < 0.01),]

library(devtools)
load_all()

infotmp <- read.table.rows(filename="/home/adai/genotyped/imputed_v3/imputed.info",
                           gsub("chr","",eQTLPUTM$snps),sep=" ")

variantsPUTM <- infotmp[as.character(gsub("chr","",eQTLPUTM$snps)),]


tmpVariansPUTM <- as.data.frame(cbind(as.character(eQTLPUTM$rs),as.character(gsub("chr","",eQTLPUTM$snps)),
                                      as.character(eQTLPUTM$gene),eQTLPUTM$pvalue))
head(tmpVariansPUTM)

tmpVariansPUTM$chr <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,2]),":")
                                    ,function(x){c(x[1])}))

tmpVariansPUTM$BP <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,2]),":")
                                   ,function(x){c(x[2])}))

tmpVariansPUTM$tmp <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,2]),":")
                                    ,function(x){c(x[3])}))

tmpVariansPUTM$Allele1 <- as.character(variantsPUTM$Al1)

tmpVariansPUTM$Allele2 <- as.character(variantsPUTM$Al2)

for (i in 1:nrow(tmpVariansPUTM))
{
  if(!is.na(tmpVariansPUTM$tmp[i]))
  {
    tmp <- unlist(strsplit(as.character(tmpVariansPUTM$tmp[i]),"_"))
    
    tmpVariansPUTM$Allele1[i] <- tmp[1]
    tmpVariansPUTM$Allele2[i] <- tmp[2]
  }
}


tmpVariansPUTM$Name <- unlist(lapply(strsplit(as.character(tmpVariansPUTM[,1]),";")
                                     ,function(x){c(x[1])}))
head(tmpVariansPUTM)

head(variantsPUTM)

tmpVariansPUTM <- tmpVariansPUTM[,c(10,5,6,8,9,3,4)]

head(tmpVariansPUTM[,c(10,5,6,8,9,3,4)])

colnames(tmpVariansPUTM) <- c("Name","CHR","BP","Allele1","Allele2","Gene","Pval")

write.csv(tmpVariansPUTM,file="eQTL.PUTM.csv",row.names=F)

## SNIG

infotmp <- read.table.rows(filename="/home/adai/genotyped/imputed_v3/imputed.info",
                           gsub("chr","",eQTLSNIG$snps),sep=" ")

variantsSNIG <- infotmp[as.character(gsub("chr","",eQTLSNIG$snps)),]


tmpVariansSNIG <- as.data.frame(cbind(as.character(eQTLSNIG$rs),as.character(gsub("chr","",eQTLSNIG$snps)),
                                      as.character(eQTLSNIG$gene),eQTLSNIG$pvalue))
head(tmpVariansSNIG)

tmpVariansSNIG$chr <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,2]),":")
                                    ,function(x){c(x[1])}))

tmpVariansSNIG$BP <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,2]),":")
                                   ,function(x){c(x[2])}))

tmpVariansSNIG$tmp <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,2]),":")
                                    ,function(x){c(x[3])}))

tmpVariansSNIG$Allele1 <- as.character(variantsSNIG$Al1)

tmpVariansSNIG$Allele2 <- as.character(variantsSNIG$Al2)

for (i in 1:nrow(tmpVariansSNIG))
{
  if(!is.na(tmpVariansSNIG$tmp[i]))
  {
    tmp <- unlist(strsplit(as.character(tmpVariansSNIG$tmp[i]),"_"))
    
    tmpVariansSNIG$Allele1[i] <- tmp[1]
    tmpVariansSNIG$Allele2[i] <- tmp[2]
  }
}


tmpVariansSNIG$Name <- unlist(lapply(strsplit(as.character(tmpVariansSNIG[,1]),";")
                                     ,function(x){c(x[1])}))
head(tmpVariansSNIG)

head(variantsSNIG)

tmpVariansSNIG <- tmpVariansSNIG[,c(10,5,6,8,9,3,4)]

colnames(tmpVariansSNIG) <- c("Name","CHR","BP","Allele1","Allele2","Gene","Pval")

write.csv(tmpVariansSNIG,file="eQTL.SNIG.csv",row.names=F)



