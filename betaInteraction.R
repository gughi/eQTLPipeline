# beta intercation test between intronic and exonic eQTLs


## calculate the beta interaction for my eQTLs
## the FUNCTION TO CALCULATE THE INTERACTION IS AT THE BOTTOM
## 27-1-14


load_all()
library(lme4)

load("data/general/overlappingGenes.rda")

## load the eQTLs
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR <=0.05),]
eQTLExonic <- eQTLPUTM[-which(eQTLPUTM$gene %in% names(listNonOve)),c("snps","gene")]
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR <=0.05),]
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
  ##snp0 <- round(snp0)
  
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



save(res,file="data/results/betaInteractionExIn5FDR.PUTM.rda")
write.csv(res,file="data/results/betaInteractionExIn5FDR.PUTM.csv")
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
  # snp0 <- round(snp0)
  
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




###############################
### TEST LINEAR MODEL
###############################


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
  ##snp0 <- round(snp0)
  
  ## markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
  
  ## rownames(sampleInfo) <- sampleInfo$A.CEL_file
  
  colnames(snp0) <- rownames(expr)
  
  f <- betaInteractionExonicIntronic(expr,snp0)
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
  
  restmp <- cbind(t(as.matrix(f)),statExonic,statIntronic)
  colnames(restmp)[1:3] <- c("beta.interaction","Std.Error.interaction","p.interaction")
  
  if (is.null(res))
  {
    res <- restmp
  } else {
    res <- rbind(res,restmp)
  }
  print(j)
  rm(statExonic,statIntronic,expr,markers,markers.info,snp0,selSamples,restmp,f)
}



save(res,file="data/results/betaInteractionExInInterTerm.SNIG.rda")
write.csv(res,file="data/results/betaInteractionExInInterTerm.SNIG.csv")
## qqplot
qq.plot(x=res$p.interaction,main="beta interaction test SNIG")



#############
### SNIG ####
#############



load("data/general/overlappingGenes.rda")

## load the eQTLs
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
eQTLExonic <- eQTLSNIG[-which(eQTLSNIG$gene %in% names(listNonOve)),c("snps","gene")]
load("data/results/finaleQTLs/intronic.Ann.SNIG.rda")
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
  ##snp0 <- round(snp0)
  
  ## markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
  
  ## rownames(sampleInfo) <- sampleInfo$A.CEL_file
  
  colnames(snp0) <- rownames(expr)
  
  f <- betaInteractionExonicIntronic(expr,snp0)
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
  
  restmp <- cbind(t(as.matrix(f)),statExonic,statIntronic)
  colnames(restmp)[1:3] <- c("beta.interaction","Std.Error.interaction","p.interaction")
  
  if (is.null(res))
  {
    res <- restmp
  } else {
    res <- rbind(res,restmp)
  }
  print(j)
  rm(statExonic,statIntronic,expr,markers,markers.info,snp0,selSamples,restmp,f)
}



save(res,file="data/results/betaInteractionExInInterTerm.SNIG.rda")
write.csv(res,file="data/results/betaInteractionExInInterTerm.SNIG.csv")
## qqplot
qq.plot(x=res$p.interaction,main="beta interaction test SNIG")










## Mike suggested to stick to the first model (mixed model) and instead of doing Q-Q plot do a volcano plot. 



load("data/results/betaInteractionExIn.PUTM.rda")
fold.threshold <- 2
p.val.threshold <- 1e-05


for(j in 1:nrow(res)){
  
  #print(expr.data.1[,gene])
  #   t.tests[[gene]] <- t.test(expr.data.1[,colnames(expr.data.1) %in% gene],expr.data.2[,colnames(expr.data.2) %in% gene])  	
  #   p.vals[[gene]] <- -log10(t.tests[[gene]]$p.value)
  #print(t.tests[[gene]]$p.value)
  beta.e <- res[j,"ge.beta"]
  beta.i <- res[j,"gi.beta"]
  res[j,"fold.change"] <- beta.i/beta.e
  if(res[j,"fold.change"] < 1){
    res[j,"fold.change"] <- -1/res[j,"fold.change"]
  }
  if(abs(res[j,"fold.change"]) > fold.threshold & res[j,"p.interaction"] < p.val.threshold){
    res[j,"colors"] <- "red"  
  }else
    res[j,"colors"] <- "black"	
  res$p.interaction

  
}

plot.title = paste0("Vulcano plot, all TFs, ",
                    tissue[1],"(left) modules vs. ",tissues[2]," tissue, fold th=",
                    fold.threshold, " and -log10(pval) th=",p.val.threshold)
pdf(paste0(plot.path,"/vulcano_all.",tissue,".pdf"),width=14,height=10)
legend.text = c("High fold and signif. p val","Nothing relevant")



par(mar=c(4, 4, 4, 4))
res <- res[-which(sign(res$ge.beta) != sign(res$gi.beta)),]
plot(res$fold.change,-log(res$p.interaction),t="p",col=res$colors,
#     main=plot.title,
     xlab="Fold change",ylab="-log10(p.val)",bg=colors,pch=21,cex=1)
,
##     cex=2*mm^2)





max(res$fold.change)


res[which(res$fold.change >100),]




if(abs(res[64,"fold.change"]) > fold.threshold & res["p.interaction"] < p.val.threshold)


tmp1 <- res[which(res[,"p.interaction"] < p.val.threshold),]




intTerm <- read.csv("data/results/betaInteractionExInInterTerm.PUTM.csv")

tmp2[-which(sign(res$ge.beta) != sign(res$gi.beta)),]

tmp2 <- intTerm[which(intTerm[,"p.interaction"] < p.val.threshold),]


length(intersect(paste0(tmp1$ge.SNP,tmp1$ge.gene),paste0(tmp2$ge.SNP,tmp2$ge.gene)))

tmp1[which(tmp1$fold.change > 2),]
  


  load("data/results/betaInteractionExIn.PUTM.rda")
  qq.plot(p.adjust(res$p.interaction,method="fdr",n=674))
         
  res$FDRInter <- p.adjust(res$p.interaction,method="fdr",n=674)
  fold.threshold <- 2
  p.val.threshold <- 0.05
  
                  
  for(j in 1:nrow(res)){
                    
    beta.e <- res[j,"ge.beta"]
    beta.i <- res[j,"gi.beta"]
    res[j,"fold.change"] <- beta.i/beta.e
    if(res[j,"fold.change"] < 1){
        res[j,"fold.change"] <- -1/res[j,"fold.change"]
    }
    if(abs(res[j,"fold.change"]) > fold.threshold & res[j,"FDRInter"] < p.val.threshold){
      res[j,"colors"] <- "red"  
    }else
      res[j,"colors"] <- "black"	
      res$p.interaction
  }

  plot.title = paste0("Vulcano plot, all TFs, ",
            tissue[1],"(left) modules vs. ",tissues[2]," tissue, fold th=",
            fold.threshold, " and -log10(pval) th=",p.val.threshold)
            pdf(paste0(plot.path,"/vulcano_all.",tissue,".pdf"),width=14,height=10)
            legend.text = c("High fold and signif. p val","Nothing relevant")
                
                  
                  
  par(mar=c(4, 4, 4, 4))
  res <- res[-which(sign(res$ge.beta) != sign(res$gi.beta)),]
  ##res <- res[- which(res$fold.change >100),]
  plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
  #     main=plot.title,
  xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR",bg=colors,pch=21,cex=1)
  
  ##     cex=2*mm^2)
                  
  dim(res[which(res[,"FDRInter"] < p.val.threshold),])
  text(x=res[which(res$FDRInter <0.05 & res$fold.change >20),"fold.change"],
       y=-log10(res[which(res$FDRInter <0.05 & res$fold.change >20),"FDRInter"]),
       labels=paste(res[which(res$FDRInter <0.05 & res$fold.change >20),c("ge.gene")]),pos=2)
  
                                
  plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
    #     main=plot.title,
    xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR",bg=colors,pch=21,cex=1)
      idx <- which(res$FDRInter <1e-15 & res$fold.change < -5)             
      text(x=res[idx,"fold.change"],
      y=-log10(res[idx,"FDRInter"]),
      labels=paste(res[idx,c("ge.gene")]),pos=2)

  
         
         plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
              #     main=plot.title,
              xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR",bg=colors,pch=21,cex=1)
         idx <- which(res$FDRInter <1e-8 & res$fold.change < -0 & res$colors =="black"  )             
         text(x=res[idx,"fold.change"],
              y=-log10(res[idx,"FDRInter"]),
              labels=paste(res[idx,c("ge.gene")]),pos=2)
         

                  





## below we try to divide the differences between exonic and intronic eQTLs
load("data/results/betaInteractionExIn.PUTM.rda")
fold.threshold <- 2
p.val.threshold <- 1e-05

for(j in 1:nrow(res)){
  
  #print(expr.data.1[,gene])
  #   t.tests[[gene]] <- t.test(expr.data.1[,colnames(expr.data.1) %in% gene],expr.data.2[,colnames(expr.data.2) %in% gene])    
  #   p.vals[[gene]] <- -log10(t.tests[[gene]]$p.value)
  #print(t.tests[[gene]]$p.value)
  beta.e <- res[j,"ge.beta"]
  beta.i <- res[j,"gi.beta"]
  res[j,"fold.change"] <- beta.i/beta.e
  if(res[j,"fold.change"] < 1){
    res[j,"fold.change"] <- -1/res[j,"fold.change"]
  }
  if(abs(res[j,"fold.change"]) > fold.threshold & res[j,"p.interaction"] < p.val.threshold){
    res[j,"colors"] <- "red"  
  }else
    res[j,"colors"] <- "black"	
  res$p.interaction
  rm(beta.e,beta.i)  
  
}


ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=res$ge.gene, mart=ensembl)

res$TSS <- sapply(res$ge.gene, function(x){getTSS(x,geneNames)})

posSNP <- unlist(lapply(strsplit(as.character(res$ge.SNP),":"),function(x){x[2]}))
DisGeneStart <- res$TSS - as.integer(posSNP)
res$DisGeneStart <- DisGeneStart
rm(DisGeneStart,posSNP)

res <- res[-which(sign(res$ge.beta) != sign(res$gi.beta)),]
plot(res$fold.change,-log(res$p.interaction),t="p",col=res$colors,
     #     main=plot.title,
     xlab="Fold change",ylab="-log10(p.val)",bg=colors,pch=21,cex=1)

## select significant and positive fold change greater than 2
idx <- which(res$p.interaction <1e-5 & res$fold.change > 2)             
positive <- res[idx,]  
res <- res[-idx,]

## select significant and positive fold change lesser than -2
idx <- which(res$p.interaction <1e-5 & res$fold.change < -2)             
negative <- res[idx,]  
res <- res[-idx,]


par(mfrow=c(1,1))
hist(negative$DisGeneStart,col='skyblue',border=F,main= "TSS Intronic vs Exonic",
     sub=paste("Intronic:",length(negative$DisGeneStart),"Exonic:",length(positive$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(negative$DisGeneStart,positive$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(positive$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(negative$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(positive$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("Exonic","Intronic"),col=c('skyblue','red'),pch=15)



par(mfrow=c(1,1))
hist(c(negative$DisGeneStart,positive$DisGeneStart),col='skyblue',border=F,main= "TSS Significant vs non-significant",
     sub=paste("Significant beta interaction:",length(c(negative$DisGeneStart,positive$DisGeneStart)),
               "non significant beta interaction:",length(res$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(c(negative$DisGeneStart,positive$DisGeneStart),res$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(c(negative$DisGeneStart,positive$DisGeneStart), adjust = 2), col = "skyblue")
lines(density(res$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("significant","non significant"),col=c('skyblue','red'),pch=15)

######################################################
## We categories the eQTLs in the beta interaction ###
######################################################

############
### PUTM ### 
############
{
load("data/results/betaInteractionExIn.PUTM.rda")
qq.plot(p.adjust(res$p.interaction,method="fdr",n=674))

res$FDRInter <- p.adjust(res$p.interaction,method="fdr",n=674)
fold.threshold <- 2
p.val.threshold <- 0.05

for(j in 1:nrow(res)){
  
   beta.e <- res[j,"ge.beta"]
   beta.i <- res[j,"gi.beta"]
   res[j,"fold.change"] <- beta.i/beta.e
   if(res[j,"fold.change"] < 1){
     res[j,"fold.change"] <- -1/res[j,"fold.change"]
   }
  
  if(res[j,"fold.change"] > fold.threshold & res[j,"FDRInter"] < p.val.threshold){
    res[j,"colors"] <- "red"
  }
  else if(res[j,"fold.change"] < -fold.threshold & res[j,"FDRInter"] < p.val.threshold){
    res[j,"colors"] <- "blue"  
  }else
    res[j,"colors"] <- "black"  
  res$p.interaction
}
par(mar=c(4, 4, 4, 4))
res <- res[-which(sign(res$ge.beta) != sign(res$gi.beta)),]
##res <- res[- which(res$fold.change >100),]
plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
     #     main=plot.title,
     xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR (PUTM)",bg=colors,pch=21,cex=1)
legend("topright",c("'positive'","'negative'","non significant"),col=c('red','blue','black'),pch=15)
##     cex=2*mm^2)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=res$ge.gene, mart=ensembl)


res$TSS <- sapply(res$ge.gene, function(x){getTSS(x,geneNames)})
posSNP <- unlist(lapply(strsplit(as.character(res$ge.SNP),":"),function(x){x[2]}))
DisGeneStart <- res$TSS - as.integer(posSNP)
res$DisGeneStart <- DisGeneStart
rm(DisGeneStart,posSNP,beta.e,beta.i,j)

## select significant and positive fold change greater than 2
idx <- which(res$FDRInter < p.val.threshold & res$fold.change > fold.threshold)             
positive <- res[idx,]  
res <- res[-idx,]

## select significant and positive fold change lesser than -2
idx <- which(res$FDRInter < p.val.threshold & res$fold.change < -fold.threshold)             
negative <- res[idx,]  
res <- res[-idx,]


par(mfrow=c(1,1))
hist(positive$DisGeneStart,col='skyblue',border=F,main= "TSS 'postitive' vs non-significant (PUTM)",
     sub=paste("Significant beta interaction:",length(positive$DisGeneStart),
               "non significant beta interaction:",length(res$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(positive$DisGeneStart,res$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(positive$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(res$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("'positive'","non significant"),col=c('skyblue','red'),pch=15)


par(mar=c(6,3,3,2))
hist(c(negative$DisGeneStart),col='skyblue',border=F,main= "TSS 'negative' vs non-significant (PUTM)",
     sub=paste("Significant beta interaction:",length(negative$DisGeneStart),
               "non significant beta interaction:",length(res$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(negative$DisGeneStart,res$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(negative$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(res$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("'negative'","non significant"),col=c('skyblue','red'),pch=15)


rm(idx,p.val.threshold)

negative$category <- "negative"
positive$category <- "positive"
res$category <- "nonSignif"

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","strand","gene_biotype"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=c(as.character(res$ge.gene),
                            as.character(positive$ge.gene),
                            as.character(negative$ge.gene)), mart=ensembl)

rownames(geneNames) <- geneNames$ensembl_gene_id

positive$geneSymbol <- geneNames[as.character(positive$ge.gene),"external_gene_id"]
positive$gene_biotype <- geneNames[as.character(positive$ge.gene),"gene_biotype"]

negative$geneSymbol <- geneNames[as.character(negative$ge.gene),"external_gene_id"]
negative$gene_biotype <- geneNames[as.character(negative$ge.gene),"gene_biotype"]

res$geneSymbol <- geneNames[as.character(res$ge.gene),"external_gene_id"]
res$gene_biotype <- geneNames[as.character(res$ge.gene),"gene_biotype"]



par(mar=c(11,5,3,2))
barplot(sort(table(negative$gene_biotype),decreasing=T)/length(negative$gene_biotype),
        col='skyblue',border=F,main= "biotype  'negative' VS 'postitive' (PUTM)",
        las=2,ylim=c(0,1))
barplot(sort(table(positive$gene_biotype),decreasing=T)/length(positive$gene_biotype),
        col=scales::alpha('red',.5),border=F,add=T,
        las=2)
legend("topright",c("'negative'","'positive'"),col=c('skyblue','red'),pch=15)

finalTable <- rbind(res[,c("ge.SNP","ge.gene", "geneSymbol","gene_biotype",
          "ge.FDR","gi.FDR","FDRInter","ge.beta","gi.beta","fold.change","category")],
      positive[,c("ge.SNP","ge.gene", "geneSymbol","gene_biotype",
             "ge.FDR","gi.FDR","FDRInter","ge.beta","gi.beta","fold.change","category")],
      negative[,c("ge.SNP","ge.gene", "geneSymbol","gene_biotype",
             "ge.FDR","gi.FDR","FDRInter","ge.beta","gi.beta","fold.change","category")])

library(doParallel)
library(foreach)

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_SNP", host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_snp")


detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(20)
clusterExport(cl, c("annSinSNP","getBM"))

registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
conse <- foreach(i=1:nrow(finalTable),.combine=rbind,.verbose=F)%dopar%annSinSNP(finalTable[i,1],ensembl)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)
colnames(conse) <- c("rsID","consequence")
finalTable$rs <- conse[,1]

write.csv(finalTable,file="data/results/finalTableBetaInteraction.PUTM.csv")
}

############
### SNIG ### 
############

{
######################################################
## We categories the eQTLs in the beta interaction ###
######################################################

load("data/results/betaInteractionExIn.SNIG.rda")

res$FDRInter <- p.adjust(res$p.interaction,method="fdr",n=398)
fold.threshold <- 2
p.val.threshold <- 0.05

for(j in 1:nrow(res)){
  
  beta.e <- res[j,"ge.beta"]
  beta.i <- res[j,"gi.beta"]
  res[j,"fold.change"] <- beta.i/beta.e
  if(res[j,"fold.change"] < 1){
    res[j,"fold.change"] <- -1/res[j,"fold.change"]
  }
  
  if(res[j,"fold.change"] > fold.threshold & res[j,"FDRInter"] < p.val.threshold){
    res[j,"colors"] <- "red"
  }
  else if(res[j,"fold.change"] < -fold.threshold & res[j,"FDRInter"] < p.val.threshold){
    res[j,"colors"] <- "blue"  
  }else
    res[j,"colors"] <- "black"  
  res$p.interaction
}




par(mar=c(4, 4, 4, 4))
res <- res[-which(sign(res$ge.beta) != sign(res$gi.beta)),]

##res <- res[- which(res$fold.change >100),]
## to visualise the volcano plot properly we remove outliers point with fold change less then -30 ( found only two points)
plot(res[-which(res$fold.change < -30),"fold.change"],
         -log10(res[-which(res$fold.change < -30),"FDRInter"]),t="p",
          col=res[-which(res$fold.change < -30),"colors"],
     #     main=plot.title,
     xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR (SNIG)",bg=colors,pch=21,cex=1)
legend("topright",c("'positive'","'negative'","non significant"),col=c('red','blue','black'),pch=15)
##     cex=2*mm^2)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=res$ge.gene, mart=ensembl)


res$TSS <- sapply(res$ge.gene, function(x){getTSS(x,geneNames)})
posSNP <- unlist(lapply(strsplit(as.character(res$ge.SNP),":"),function(x){x[2]}))
DisGeneStart <- res$TSS - as.integer(posSNP)
res$DisGeneStart <- DisGeneStart
rm(DisGeneStart,posSNP,beta.e,beta.i,j)

## select significant and positive fold change greater than 2
idx <- which(res$FDRInter < p.val.threshold & res$fold.change > fold.threshold)             
positive <- res[idx,]  
res <- res[-idx,]

## select significant and positive fold change lesser than -2
idx <- which(res$FDRInter < p.val.threshold & res$fold.change < -fold.threshold)             
negative <- res[idx,]  
res <- res[-idx,]


par(mfrow=c(1,1))
hist(positive$DisGeneStart,col='skyblue',border=F,main= "TSS 'postitive' vs non-significant (SNIG)",
     sub=paste("Significant beta interaction:",length(positive$DisGeneStart),
               "non significant beta interaction:",length(res$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(positive$DisGeneStart,res$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(positive$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(res$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("'positive'","non significant"),col=c('skyblue','red'),pch=15)

par(mar=c(6,3,3,2))
hist(c(negative$DisGeneStart),col='skyblue',border=F,main= "TSS 'negative' vs non-significant (SNIG)",
     sub=paste("Significant beta interaction:",length(negative$DisGeneStart),
               "non significant beta interaction:",length(res$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(negative$DisGeneStart,res$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(negative$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(res$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("'negative'","non significant"),col=c('skyblue','red'),pch=15)

rm(idx,p.val.threshold)

negative$category <- "negative"
positive$category <- "positive"
res$category <- "nonSignif"

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","strand","gene_biotype"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=c(as.character(res$ge.gene),
                            as.character(positive$ge.gene),
                            as.character(negative$ge.gene)), mart=ensembl)

rownames(geneNames) <- geneNames$ensembl_gene_id

positive$geneSymbol <- geneNames[as.character(positive$ge.gene),"external_gene_id"]
positive$gene_biotype <- geneNames[as.character(positive$ge.gene),"gene_biotype"]

negative$geneSymbol <- geneNames[as.character(negative$ge.gene),"external_gene_id"]
negative$gene_biotype <- geneNames[as.character(negative$ge.gene),"gene_biotype"]

res$geneSymbol <- geneNames[as.character(res$ge.gene),"external_gene_id"]
res$gene_biotype <- geneNames[as.character(res$ge.gene),"gene_biotype"]



par(mar=c(11,5,3,2))
barplot(sort(table(negative$gene_biotype),decreasing=T)/length(negative$gene_biotype),
        col='skyblue',border=F,main= "biotype  'negative' VS 'postitive' (SNIG)",
        las=2,ylim=c(0,1),ylab=)

tmpTab <- table(positive$gene_biotype)
tmpTab <- c(tmpTab["protein_coding"],
            lincRNA=0,
            tmpTab["pseudogene"])
barplot(tmpTab/length(positive$gene_biotype),
        col=scales::alpha('red',.5),border=F,add=T,
        las=2)
legend("topright",c("'negative'","'positive'"),col=c('skyblue','red'),pch=15)

finalTable <- rbind(res[,c("ge.SNP","ge.gene", "geneSymbol","gene_biotype",
                           "ge.FDR","gi.FDR","FDRInter","ge.beta","gi.beta","fold.change","category")],
                    positive[,c("ge.SNP","ge.gene", "geneSymbol","gene_biotype",
                                "ge.FDR","gi.FDR","FDRInter","ge.beta","gi.beta","fold.change","category")],
                    negative[,c("ge.SNP","ge.gene", "geneSymbol","gene_biotype",
                                "ge.FDR","gi.FDR","FDRInter","ge.beta","gi.beta","fold.change","category")])

library(doParallel)
library(foreach)

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_SNP", host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_snp")


detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(20)
clusterExport(cl, c("annSinSNP","getBM"))

registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
conse <- foreach(i=1:nrow(finalTable),.combine=rbind,.verbose=F)%dopar%annSinSNP(finalTable[i,1],ensembl)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)
colnames(conse) <- c("rsID","consequence")
finalTable$rs <- conse[,1]

write.csv(finalTable,file="data/results/finalTableBetaInteraction.SNIG.csv")

}






