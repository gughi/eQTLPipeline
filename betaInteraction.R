# beta intercation test between intronic and exonic eQTLs


## calculate the beta interaction for my eQTLs
## the FUNCTION TO CALCULATE THE INTERACTION IS AT THE BOTTOM
## 27-1-14


load_all()
library(lme4)

load("data/general/overlappingGenes.rda")

## load the eQTLs
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLExonic <- eQTLPUTM[-which(eQTLPUTM$gene %in% names(listNonOve)),c("snps","gene")]
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLIntronic <- eQTLPUTM[-which(eQTLPUTM$gene %in% names(listNonOve)),c("snps","gene")]
rm(eQTLPUTM)

## get the beta interaction for intronic and exonic
eQTLExonic$type <- "exonic"
eQTLIntronic$type <- "intronic"

eQTLs <- rbind(eQTLExonic,eQTLIntronic)
rm(eQTLExonic,eQTLIntronic)


## load the expression(residual corrected values)
load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")
exprExonic <- resids
load("data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda")
exprIntronic <- resids
rm(resids,listNonOve)

eQTLs <- eQTLs[which(as.character(eQTLs$gene) %in% as.character(colnames(exprIntronic))),]
eQTLs <- eQTLs[which(as.character(eQTLs$gene) %in% as.character(colnames(exprExonic))),]

## remove eQTLs that are duplicated
eQTLs <- eQTLs[-which(duplicated(eQTLs[,1:2])),]

res <- NULL
for (j in 1:nrow(eQTLs))
{
  

  ## load the expression for each individual gene for exonic and intronic
  exprE <- exprExonic[,as.character(eQTLs$gene[j])]
  exprI <- exprIntronic[,as.character(eQTLs$gene[j])]
  
  ## add to the expression the tissue information or the exonic/intronic information
  expr <- cbind(exprE,"Exonic")
  exprtmp <- cbind(exprI,"Intronic")
  expr <- rbind(expr,exprtmp)
  rm(exprtmp,exprE,exprI)
  
  ## load the snps
  load(paste0("/home/seb/eQTL/snps/byGene/",as.character(eQTLs$gene[j]),".rda"))
  
  ## load the samples information
  load("data/general/sampleInfo.rda")
  
  
  ## change the samples names and remove outlier
  sampleInfo$U.SD_No <- gsub("/","_",sampleInfo$U.SD_No)
#   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "95_40"),]
#   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "09_50"),]
#   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "034_07"),]
#   
  rownames(sampleInfo) <- sampleInfo$A.CEL_file
  selSamples <- sampleInfo[as.character(rownames(expr)),1]
  
  ##snp0 <- markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
  
  snp0 <- markers[as.character(eQTLs$snps[j]),selSamples]
  snp0 <- round(snp0)
  
  ## markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
  
  ## rownames(sampleInfo) <- sampleInfo$A.CEL_file
  
  colnames(snp0) <- rownames(expr)
    
  f <- betaInteraction(expr,snp0)
  ## collect summary statisctics from gene exonic
  cmd <- paste0("grep -w '",as.character(eQTLs$snps[j]),"' data/results/genic/geneExons/fullResults/PUTM/",as.character(eQTLs$gene[j]))
  statExonic <- read.delim(pipe(cmd),header=F)
  rm(cmd)
  colnames(statExonic) <- c("ge.SNP","ge.gene","ge.beta","ge.t-stat","ge.p-value","ge.FDR") 
  ## collect summary statisctics from gene intronic
  cmd <- paste0("grep -w '",as.character(eQTLs$snps[j]),"' data/results/genic/geneIntronic/fullResults/PUTM/",as.character(eQTLs$gene[j]))
  statIntronic <- read.delim(pipe(cmd),header=F)
  rm(cmd)
  colnames(statIntronic) <- c("gi.SNP","gi.gene","gi.beta","gi.t-stat","gi.p-value","gi.FDR") 
  
  restmp <- cbind(f,statExonic,statIntronic)
  colnames(restmp)[1] <- "p.interaction"
  
  if (is.null(res))
  {
    res <- restmp
  } else {
    res <- rbind(res,restmp)
  }
  rm(statExonic,statIntronic,expr,markers,markers.info,snp0,selSamples,restmp,f)
}

save(res,file="data/results/betaInteractionExIn.PUTM.rda")
write.csv(res,file="data/results/betaInteractionExIn.PUTM.csv")
## qqplot
qq.plot(x=res$p.interaction,main="beta interaction test PUTM")



load_all()
library(lme4)

load("data/general/overlappingGenes.rda")

## load the eQTLs
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
eQTLExonic <- eQTLSNIG[-which(eQTLSNIG$gene %in% names(listNonOve)),c("snps","gene")]
eQTLSNIG <- read.csv("data/results/finaleQTLs/intronic.Ann.SNIG.csv")
eQTLIntronic <- eQTLSNIG[-which(eQTLSNIG$gene %in% names(listNonOve)),c("snps","gene")]
rm(eQTLSNIG)

## get the beta interaction for intronic and exonic
eQTLExonic$type <- "exonic"
eQTLIntronic$type <- "intronic"

eQTLs <- rbind(eQTLExonic,eQTLIntronic)
rm(eQTLExonic,eQTLIntronic)


## load the expression(residual corrected values)
load("data/expr/normalisedCounts/genic/geneExons/resids.SNIG.rda")
exprExonic <- resids
load("data/expr/normalisedCounts/genic/geneIntronic/resids.SNIG.rda")
exprIntronic <- resids
rm(resids,listNonOve)

eQTLs <- eQTLs[which(as.character(eQTLs$gene) %in% as.character(colnames(exprIntronic))),]
eQTLs <- eQTLs[which(as.character(eQTLs$gene) %in% as.character(colnames(exprExonic))),]

## remove eQTLs that are duplicated
eQTLs <- eQTLs[-which(duplicated(eQTLs[,1:2])),]

res <- NULL
for (j in 1:nrow(eQTLs))
{
  
  
  ## load the expression for each individual gene for exonic and intronic
  exprE <- exprExonic[,as.character(eQTLs$gene[j])]
  exprI <- exprIntronic[,as.character(eQTLs$gene[j])]
  
  ## add to the expression the tissue information or the exonic/intronic information
  expr <- cbind(exprE,"Exonic")
  exprtmp <- cbind(exprI,"Intronic")
  expr <- rbind(expr,exprtmp)
  rm(exprtmp,exprE,exprI)
  
  ## load the snps
  load(paste0("/home/seb/eQTL/snps/byGene/",as.character(eQTLs$gene[j]),".rda"))
  
  ## load the samples information
  load("data/general/sampleInfo.rda")
  
  
  ## change the samples names and remove outlier
  sampleInfo$U.SD_No <- gsub("/","_",sampleInfo$U.SD_No)
  #   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "95_40"),]
  #   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "09_50"),]
  #   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "034_07"),]
  #   
  rownames(sampleInfo) <- sampleInfo$A.CEL_file
  selSamples <- sampleInfo[as.character(rownames(expr)),1]
  
  ##snp0 <- markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
  
  snp0 <- markers[as.character(eQTLs$snps[j]),selSamples]
  snp0 <- round(snp0)
  
  ## markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
  
  ## rownames(sampleInfo) <- sampleInfo$A.CEL_file
  
  colnames(snp0) <- rownames(expr)
  
  f <- betaInteraction(expr,snp0)
  ## collect summary statisctics from gene exonic
  cmd <- paste0("grep -w '",as.character(eQTLs$snps[j]),"' data/results/genic/geneExons/fullResults/SNIG/",as.character(eQTLs$gene[j]))
  statExonic <- read.delim(pipe(cmd),header=F)
  rm(cmd)
  colnames(statExonic) <- c("ge.SNP","ge.gene","ge.beta","ge.t-stat","ge.p-value","ge.FDR") 
  ## collect summary statisctics from gene intronic
  cmd <- paste0("grep -w '",as.character(eQTLs$snps[j]),"' data/results/genic/geneIntronic/fullResults/SNIG/",as.character(eQTLs$gene[j]))
  statIntronic <- read.delim(pipe(cmd),header=F)
  rm(cmd)
  colnames(statIntronic) <- c("gi.SNP","gi.gene","gi.beta","gi.t-stat","gi.p-value","gi.FDR") 
  
  restmp <- cbind(f,statExonic,statIntronic)
  colnames(restmp)[1] <- "p.interaction"
  
  if (is.null(res))
  {
    res <- restmp
  } else {
    res <- rbind(res,restmp)
  }
  rm(statExonic,statIntronic,expr,markers,markers.info,snp0,selSamples,restmp,f)
}

save(res,file="data/results/betaInteractionExIn.SNIG.rda")
write.csv(res,file="data/results/betaInteractionExIn.SNIG.csv")
## qqplot
qq.plot(x=res$p.interaction,main="beta interaction test SNIG")



















