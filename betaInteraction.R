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
load("data/results/betaInteraction/betaInteractionExIn.PUTM.rda")
qq.plot(p.adjust(res$p.interaction,method="fdr",n=674),main="beta interaction test FDR values(PUTM)")
qq.plot(res$p.interaction,main="beta interaction test p-values(PUTM)")

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

par(mar=c(6,3,3,2))
hist(c(negative$DisGeneStart),col='skyblue',border=F,main= "TSS  'negative' vs 'postive' (PUTM)",
     sub=paste("negative beta interaction:",length(negative$DisGeneStart),
               "positive beta interaction:",length(positive$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(negative$DisGeneStart,positive$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(positive$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
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

## compare the three different categories.
finalTab <- read.csv("data/results/finalTableBetaInteraction.PUTM.csv")

par(mar=c(11,5,3,2))
nonsig <- finalTab[which(finalTab$category %in% "nonSignif"),
                   c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
nonsig <- unique(nonsig[,2:5])
positive <- finalTab[which(finalTab$category %in% "positive"),
                   c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
positive <- unique(positive[,2:5])
negative <- finalTab[which(finalTab$category %in% "negative"),
                     c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
negative <- unique(negative[,2:5])

counts <- rbind(nonsig=sort(table(nonsig$gene_biotype),decreasing=T)/length(nonsig$gene_biotype),
                positive=sort(table(positive$gene_biotype),decreasing=T)/length(positive$gene_biotype),
                negative=sort(table(negative$gene_biotype),decreasing=T)/length(negative$gene_biotype))

barplot(counts, main="Biotype PUTM",
        col=1:3,
        legend = rownames(counts), beside=TRUE,las=2,ylim=c(0,1))

}

############
### SNIG ### 
############

{
######################################################
## We categories the eQTLs in the beta interaction ###
######################################################

load("data/results/betaInteraction/betaInteractionExIn.SNIG.rda")
qq.plot(p.adjust(res$p.interaction,method="fdr",n=674),main="beta interaction test FDR values(SNIG)")
qq.plot(res$p.interaction,main="beta interaction test p-values(SNIG)")


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


finalTab <- read.csv("data/results/finalTableBetaInteraction.SNIG.csv")

par(mar=c(11,5,3,2))
nonsig <- finalTab[which(finalTab$category %in% "nonSignif"),
                   c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
nonsig <- unique(nonsig[,2:5])
positive <- finalTab[which(finalTab$category %in% "positive"),
                     c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
positive <- unique(positive[,2:5])
negative <- finalTab[which(finalTab$category %in% "negative"),
                     c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
negative <- unique(negative[,2:5])

counts <- rbind(nonsig=sort(table(nonsig$gene_biotype),decreasing=T)/length(nonsig$gene_biotype),
                positive=sort(table(positive$gene_biotype),decreasing=T)/length(positive$gene_biotype),
                negative=sort(table(negative$gene_biotype),decreasing=T)/length(negative$gene_biotype))

barplot(counts, main="Biotype SNIG",
        col=1:3,
        legend = rownames(counts), beside=TRUE,las=2,ylim=c(0,1))



}

###################
### SNIG + PUTM ### 
###################

{
load("data/results/betaInteraction/betaInteractionExIn.PUTM.rda")
tmp <-res 
load("data/results/betaInteraction/betaInteractionExIn.SNIG.rda")
res <- rbind(res,tmp)
rm(tmp)
qq.plot(p.adjust(res$p.interaction,method="fdr",n=nrow(res)),main="beta interaction test FDR values(SNIG+PUTM)")
qq.plot(res$p.interaction,main="beta interaction test p-values(SNIG+PUTM)")

res$FDRInter <- p.adjust(res$p.interaction,method="fdr",n=nrow(res))
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
res <- res[-which(res$fold.change < -40),]
plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
     #     main=plot.title,
     xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR (PUTM+SNIG)",bg=colors,pch=21,cex=1)
legend("topright",c("'positive'","'negative'","non significant"),col=c('red','blue','black'),pch=15)
##     cex=2*mm^2)
plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
     #     main=plot.title,
     xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR",bg=colors,pch=21,cex=1)

## plot an example of negative 
idx <- which(res$FDRInter <1e-15 & res$fold.change < -5)             
text(x=res[idx,"fold.change"],
     y=-log10(res[idx,"FDRInter"]),
     labels=paste(res[idx,c("ge.gene")]),pos=2,col="blue")

plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
     #     main=plot.title,
     xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR",bg=colors,pch=21,cex=1)
## positive example
text(x=res[which(res$FDRInter <0.05 & res$fold.change >20),"fold.change"],
     y=-log10(res[which(res$FDRInter <0.05 & res$fold.change >20),"FDRInter"]),
     labels=paste(res[which(res$FDRInter <0.05 & res$fold.change >20),c("ge.gene")]),pos=2,col="red")


plot(res$fold.change,-log10(res$FDRInter),t="p",col=res$colors,
     #     main=plot.title,
     xlab="Fold change",ylab="-log10(FDR)",main="Volcano plot FDR",bg=colors,pch=21,cex=1)
idx <- which(res$FDRInter <1e-8 & res$fold.change < -0 & res$colors =="black"  )             
text(x=res[idx,"fold.change"],
     y=-log10(res[idx,"FDRInter"]),
     labels=paste(res[idx,c("ge.gene")]),pos=2,col="black")

rm(beta.e,beta.i,fold.threshold,j,p.val.threshold,idx)




res <- read.csv("data/results/finalTableBetaInteraction.SNIG.csv")

## get TSS and TES
# res <- read.csv("data/results/finalTableBetaInteraction.PUTM.csv")
# tmp <-res 
# res <- read.csv("data/results/finalTableBetaInteraction.SNIG.csv")
# res <- rbind(res,tmp)
# rm(tmp)


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

nonsig <- res[which(res$category %in% "nonSignif"),
                   c("ge.SNP","ge.gene","gene_biotype","DisGeneStart")]
positive <- res[which(res$category %in% "positive"),
                     c("ge.SNP","ge.gene","gene_biotype","DisGeneStart")]
negative <- res[which(res$category %in% "negative"),
                     c("ge.SNP","ge.gene","gene_biotype","DisGeneStart")]

par(mar=c(6,5,3,2))
hist(c(negative$DisGeneStart),col='skyblue',border=F,main= "Distance TSS 'negative' vs non-significant (PUTM+SNIG)",
     sub=paste("Significant beta interaction:",length(negative$DisGeneStart),
               "non significant beta interaction:",length(res$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(negative$DisGeneStart,res$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(negative$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(res$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("'negative'","non significant"),col=c('skyblue','red'),pch=15)

hist(c(positive$DisGeneStart),col='skyblue',border=F,main= "Distance TSS 'positive' vs non-significant (PUTM+SNIG)",
     sub=paste("Significant beta interaction:",length(positive$DisGeneStart),
               "non significant beta interaction:",length(res$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(positive$DisGeneStart,res$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(positive$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(res$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("'positive'","non significant"),col=c('skyblue','red'),pch=15)

rm(negative,positive,nonsig)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","strand"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=res$ge.gene, mart=ensembl)

res$TES <- sapply(res$ge.gene, function(x){getTES(x,geneNames)})
head(res)

posSNP <- unlist(lapply(strsplit(as.character(res$ge.SNP),":"),function(x){x[2]}))
DisGeneEnd <- res$TES - as.integer(posSNP)
res$DisGeneEnd <- DisGeneEnd
rm(DisGeneEnd,posSNP)

nonsig <- res[which(res$category %in% "nonSignif"),
              c("ge.SNP","ge.gene","gene_biotype","DisGeneEnd")]
positive <- res[which(res$category %in% "positive"),
                c("ge.SNP","ge.gene","gene_biotype","DisGeneEnd")]
negative <- res[which(res$category %in% "negative"),
                c("ge.SNP","ge.gene","gene_biotype","DisGeneEnd")]

par(mar=c(6,5,3,2))
hist(c(negative$DisGeneEnd),col='skyblue',border=F,main= "Distance TES 'negative' vs non-significant (PUTM+SNIG)",
     sub=paste("Significant beta interaction:",length(negative$DisGeneEnd),
               "non significant beta interaction:",length(res$DisGeneEnd)),
     xlab=paste("KS pvalue:",ks.test(negative$DisGeneEnd,res$DisGeneEnd)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneEnd,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(negative$DisGeneEnd, adjust = 2), col = "skyblue")
lines(density(res$DisGeneEnd, adjust = 2),col = "red")
legend("topright",c("'negative'","non significant"),col=c('skyblue','red'),pch=15)

hist(c(positive$DisGeneEnd),col='skyblue',border=F,main= "Distance TES 'positive' vs non-significant (PUTM+SNIG)",
     sub=paste("Significant beta interaction:",length(positive$DisGeneEnd),
               "non significant beta interaction:",length(res$DisGeneEnd)),
     xlab=paste("KS pvalue:",ks.test(positive$DisGeneEnd,res$DisGeneEnd)$p.value),freq=FALSE,breaks = 40)
hist(res$DisGeneEnd,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(positive$DisGeneEnd, adjust = 2), col = "skyblue",freq=FALSE)
lines(density(res$DisGeneEnd, adjust = 2), col = "red")
legend("topright",c("'positive'","non significant"),col=c('skyblue','red'),pch=15)

rm(geneNames,negative,positive,nonsig,DisGeneEnd,ensembl)
#### biotype

par(mar=c(11,5,3,2))
nonsig <- res[which(as.character(res$category) %in% "nonSignif"),
                   c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
nonsig <- unique(nonsig[,2:5])
positive <- res[which(as.character(res$category) %in% "positive"),
                     c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
positive <- unique(positive[,2:5])
negative <- res[which(as.character(res$category) %in% "negative"),
                     c("ge.SNP","ge.gene","gene_biotype","ge.FDR","gi.FDR")]
negative <- unique(negative[,2:5])



counts <- rbind(nonsig=sort(table(nonsig$gene_biotype),decreasing=T)/length(nonsig$gene_biotype),
                positive=sort(table(positive$gene_biotype),decreasing=T)/length(positive$gene_biotype),
                negative=sort(table(negative$gene_biotype),decreasing=T)/length(negative$gene_biotype))

counts
table(negative$gene_biotype)
barplot(counts, main="Biotype PUTM+SNIG",
        col=1:3,
        legend = rownames(counts), beside=TRUE,las=2,ylim=c(0,1))

variantAnnoNega <- read.delim("data/results/VEP/PUTM_negative.txt")
variantAnnoPos <- read.delim("data/results/VEP/PUTM_positive.txt")
variantAnnoNonSig <- read.delim("data/results/VEP/PUTM_nonSignificant.txt")

variantAnnoNega <- rbind(variantAnnoNega,read.delim("data/results/VEP/SNIG_negative.txt"))
variantAnnoPos <- rbind(variantAnnoPos,read.delim("data/results/VEP/SNIG_positive.txt"))
variantAnnoNonSig <- rbind(variantAnnoNonSig,read.delim("data/results/VEP/SNIG_nonSignificant.txt"))


consNeg <- unlist(lapply(strsplit(as.character(variantAnnoNega$Consequence),","),function(x){x[1]}))
consPos <- unlist(lapply(strsplit(as.character(variantAnnoPos$Consequence),","),function(x){x[1]}))
consNonSig <- unlist(lapply(strsplit(as.character(variantAnnoNonSig$Consequence),","),function(x){x[1]}))

rm(variantAnnoNega,variantAnnoNonSig,variantAnnoPos)
nam <- names(sort(table(consNonSig),decreasing=T)/length(consNonSig))

nonsig <- sort(table(consNonSig),decreasing=T)/length(consNonSig)
positive <- sort(table(consPos),decreasing=T)/length(consPos)
negative <- table(consNeg)/length(consNeg)


counts <- rbind(nonsig=nonsig[nam],
                positive=positive[nam],
                negative=negative[nam])


par(mar=c(15,5,3,2))
barplot(counts, main="Variant Consequences",
        col=1:3,
        legend = rownames(counts), beside=TRUE,las=2,ylim=c(0,0.6))


}

## 26/1/2016

## follow Mike's suggestions: do the vulcano plot but instead of doing it with fold change do it with the betas in this way:

##1) Plot using a delta as follows (rather than fold-change):
##  If beta(exonic) is positive:
##  Delta = beta(exonic) - beta(intronic)
## If beta(exonic) is negative:
##   Delta = beta(intronic) - beta(exonic)

## This will mean that point on the right-hand-side of the plot are when beta(exonic) is "more extreme" (further from zero) than beta(intronic), and conversely for left-hand-side.


load("data/results/betaInteraction/betaInteractionExIn.PUTM.rda")
tmp <-res 
load("data/results/betaInteraction/betaInteractionExIn.SNIG.rda")
res <- rbind(res,tmp)
rm(tmp)

res$FDRInter <- p.adjust(res$p.interaction,method="fdr",n=nrow(res))
## calculated delta 

for(j in 1:nrow(res)){
  

  if(res[j,"ge.beta"] >= 0){
    res[j,"delta"] <- res[j,"ge.beta"] - res[j,"gi.beta"]
  }
  else
  {
    res[j,"delta"] <- res[j,"gi.beta"] - res[j,"ge.beta"]  
  }
}
par(mar=c(4, 4, 4, 4))
plot(res$delta,-log10(res$FDRInter),t="p",
     #     main=plot.title,
     xlab="delta",ylab="-log10(FDR)",main="Volcano plot FDR (PUTM+SNIG)",pch=21,cex=1)


#     cex=2*mm^2)

## For windows
## install.packages("R.utils")
library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)$pathname
setwd(dirname(path))
rm(path)

variantAnnoNega <- read.delim("data/results/VEP/PUTM_negative.txt")
variantAnnoPos <- read.delim("data/results/VEP/PUTM_positive.txt")
variantAnnoNonSig <- read.delim("data/results/VEP/PUTM_nonSignificant.txt")

variantAnnoNega <- rbind(variantAnnoNega,read.delim("data/results/VEP/SNIG_negative.txt"))
variantAnnoPos <- rbind(variantAnnoPos,read.delim("data/results/VEP/SNIG_positive.txt"))
variantAnnoNonSig <- rbind(variantAnnoNonSig,read.delim("data/results/VEP/SNIG_nonSignificant.txt"))


consNeg <- unlist(lapply(strsplit(as.character(variantAnnoNega$Consequence),","),function(x){x[1]}))
consPos <- unlist(lapply(strsplit(as.character(variantAnnoPos$Consequence),","),function(x){x[1]}))
consNonSig <- unlist(lapply(strsplit(as.character(variantAnnoNonSig$Consequence),","),function(x){x[1]}))

rm(variantAnnoNega,variantAnnoNonSig,variantAnnoPos)
nam <- names(sort(table(consNonSig),decreasing=T)/length(consNonSig))

nonsig <- sort(table(consNonSig),decreasing=T)
positive <- sort(table(consPos),decreasing=T)
negative <- table(consNeg)


counts <- rbind(nonsig=nonsig[nam],
                positive=positive[nam],
                negative=negative[nam])


## Percentages
# nonsig <- sort(table(consNonSig),decreasing=T)/length(consNonSig)
# positive <- sort(table(consPos),decreasing=T)/length(consPos)
# negative <- table(consNeg)/length(consNeg)

# counts <- rbind(nonsig=nonsig[nam],
#                 positive=positive[nam],
#                 negative=negative[nam])

counts[is.na(counts)]=0
counts <- rbind(counts,total=apply(counts,2,sum))
counts <- cbind(counts,total=apply(counts,1,sum))
ftable(t(counts))
##counts <- t(counts) 


## we test the groups we have defined as non-sig,positive and negative
pval <- NULL
for(i in 1:(nrow(counts)-1)){
  
  pval[i] <- chisq.test(counts[4,1:ncol(counts)-1],counts[i,1:ncol(counts)-1])$p.value
  
}

rm(pval,i)

## with intronic as reference
pval <- NULL
for(i in 2:(ncol(counts)-1)){
  
  pval[i] <- chisq.test(counts[1:nrow(counts)-1,1],counts[1:(nrow(counts)-1),i])$p.value
  
}

rm(pval,i)

pval <- NULL
for(i in 1:(ncol(counts)-1)){
  
  pval[i] <- chisq.test(counts[1:nrow(counts)-1,13],counts[1:(nrow(counts)-1),i])$p.value
  
}



  

## re-do contingency table using the deltas to divide the groups instead of using the 
## fold-change

{

  library(R.utils)
  path <- readWindowsShortcut("data.lnk", verbose=FALSE)
  setwd(dirname(path$networkPathname))
  rm(path)
  
  load("data/results/betaInteraction/betaInteractionExIn.PUTM.rda")
  tmp <-res 
  load("data/results/betaInteraction/betaInteractionExIn.SNIG.rda")
  res <- rbind(res,tmp)
  rm(tmp)
  
  res$FDRInter <- p.adjust(res$p.interaction,method="fdr",n=nrow(res))
  
  ## adding the one condition more for when the betas have oppositive signs
  
  for(j in 1:nrow(res)){
    
    if(sign(res[j,"ge.beta"])!=sign(res[j,"gi.beta"])){
      res[j,"delta"] <- abs(res[j,"ge.beta"]) - abs(res[j,"gi.beta"])
    }else if(res[j,"ge.beta"] >= 0){
      res[j,"delta"] <- res[j,"ge.beta"] - res[j,"gi.beta"]
    }else
    {
      res[j,"delta"] <- res[j,"gi.beta"] - res[j,"ge.beta"]  
    }
  }
  
  par(mar=c(4, 4, 4, 4))
  plot(res$delta,-log10(res$FDRInter),t="p",
       #     main=plot.title,
       xlab="delta",ylab="-log10(FDR)",main="Volcano plot FDR (PUTM+SNIG)",pch=21,cex=1)

  ##
  fdr.threshold <- 0.05
  
  for(j in 1:nrow(res)){
    
    if(res[j,"delta"] >= 0 & res[j,"FDRInter"] < fdr.threshold){
      res[j,"colors"] <- "red"
    }
    else if(res[j,"delta"] < 0 & res[j,"FDRInter"] < fdr.threshold){
      res[j,"colors"] <- "blue"  
    }else
      res[j,"colors"] <- "black"  
    res$p.interaction
  }
  
  plot(res$delta,-log10(res$FDRInter),t="p",col=res$colors,
       #     main=plot.title,
       xlab="Delta",ylab="-log10(FDR)",main="Volcano plot FDR (PUTM+SNIG)",bg=colors,pch=21,cex=1)
  legend("topright",c("'positive'","'negative'","non significant"),col=c('red','blue','black'),pch=15)
  rm(ensembl,j)
  
  library(biomaRt)
  ## annotation of biotype
  ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                     dataset="hsapiens_gene_ensembl")
  
  geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","strand","gene_biotype"),
                     verbose = T,
                     filters="ensembl_gene_id",
                     values=c(as.character(res$ge.gene)), mart=ensembl)
  
  
  
  rm(ensembl)
  
  rownames(geneNames) <- geneNames$ensembl_gene_id
  
  res <- cbind(res,geneNames[as.character(res$ge.gene),c("external_gene_id","gene_biotype")])
  
## to load the libraries
  library(devtools)
  setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
  load_all()
  path <- readWindowsShortcut("data.lnk", verbose=FALSE)
  setwd(dirname(path$networkPathname))
  rm(path)
  
  ## get the right TES look at the strand too
  res$TES <- sapply(res$ge.gene, function(x){getTES(x,geneNames)})

  ## get the distance position to the TES
  posSNP <- unlist(lapply(strsplit(as.character(res$ge.SNP),":"),function(x){x[2]}))
  DisGeneEnd <- res$TES - as.integer(posSNP)
  res$DisGeneEnd <- DisGeneEnd
  rm(DisGeneEnd,posSNP)
  
  ## we now get the 
  load("data/general/chrpos2rsid_UKBEC.rda")
  save(chrpos2rsid,file="data/general/chrpos2rsid_UKBEC.rda")
  
  tmp<- sapply(paste0("^",gsub("chr","",as.character(res$ge.SNP)),"$"),function(y) grep(y,x=as.character(chrpos2rsid$SNP)))
  
  tmp <- chrpos2rsid[tmp,]
  res$rsID <- tmp$rsid  
  
  save(res,file="data/results/betaInteraction/betaInteractionExIn.rda")
  load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

  ## We obtain the rs number for indels 
  
  library(doParallel)
  library(foreach)
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
  conse <- foreach(i=1:length(res[which(is.na(res$rsID)),"ge.SNP"]),.combine=rbind,.verbose=F)%dopar%annSinSNP(res[which(is.na(res$rsID)),"ge.SNP"][i],ensembl)
  ##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
  end <- Sys.time()
  end-start
  stopCluster(cl)

  conse <- cbind(as.character(res[which(is.na(res$rsID)),"ge.SNP"]),conse)
  
  res$rsID <- as.character(res$rsID)
  
  res[which(as.character(res$ge.SNP) %in% as.character(conse[,1])),"rsID"] <- as.character(conse[,2])
  
  head(res[which(as.character(res$ge.SNP) %in% as.character(conse[,1])),])
  
  
  save(res,file="data/results/betaInteraction/betaInteractionExIn.rda")
  
  ## we select now the different categories exonic intronic
  
  load(file="data/results/betaInteraction/betaInteractionExIn.rda")
  
  ## get the information the VEP information
  
  ## first we need to get the input file for VEP
  head(res)
  
  load("data/general/imputed.info.rda")
  
  tabVariants <- imputed.info[which(imputed.info$marker %in% as.character(gsub("chr","",res$ge.SNP))),]
  
  ## select the indels
  tabVariants$type <- unlist(lapply(strsplit(as.character(tabVariants$marker),":"),function(x){length(x)})) 
  
  tabVariants$type[which(tabVariants$type==2)] = "SNP"
  tabVariants$type[which(tabVariants$type==3)] = "indel"

    
  tabVariants <- cbind(unlist(lapply(strsplit(gsub("chr","",tabVariants$marker),":"),function(x){x[1]})),
                       unlist(lapply(strsplit(gsub("chr","",tabVariants$marker),":"),function(x){x[2]})),
                       tabVariants)
    
  tabVariants <- tabVariants[,1:5]
  
  colnames(tabVariants) <- c("#CHROM","POS","ID","REF","ALT")

  write("##fileformat=VCFv4.0",file="data/results/VEP/toAnnotate.vcf")
  write.table(tabVariants,file="data/results/VEP/toAnnotate.vcf",row.names=F,append = T,quote = F)
  ##rm(tabVariants)
  

  tabVariants <- tabVariants[which(tabVariants$type == "indel"),]
  tabVariants$Al1 <- as.character(tabVariants$Al1)
  tabVariants$Al1[which(as.character(tabVariants$Al1)=="R")]="."
  tabVariants$Al2 <- as.character(tabVariants$Al2)
  tabVariants$Al2[which(as.character(tabVariants$Al2)=="R")]="."
  
  apply(tabVariants,1,function(x){
    
    tabVariants$Al1[]
    
  })
  
  
  
  
  
  
  
  ## we annotated teh variants using caprica. in windows it wasn't possible
  ## teh command used is the following 
  
  ## perl /home/seb/softwares/VEP/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl -i tab.vcf -o testRs.txt.out --offline --dir_cache "/home/seb/.vep/" --assembly "GrCh37" --cache_version "72" --species "homo_sapiens" --symbol --sift b --poly b --hgvs --regulatory --biotype --gmaf --maf_1kg --maf_esp --pick --pick_order rank,canonical,tsl,biotype,length --fork 4
  
  
  
  
  
  
  
  write.table(unlist(lapply(strsplit(as.character(res[which(as.character(res$colors) %in% "red"),"rsID"]),";"),function(x){x})),
              file="data/results/VEP/positive.delta.txt",col.names=F,row.names=F)
  write.table(unlist(lapply(strsplit(as.character(res[which(as.character(res$colors) %in% "blue"),"rsID"]),";"),function(x){x})),
              file="data/results/VEP/negative.delta.txt",col.names=F,row.names=F)
  write.table(unlist(lapply(strsplit(as.character(res[which(as.character(res$colors) %in% "black"),"rsID"]),";"),function(x){x})),
              file="data/results/VEP/nonSig.delta.txt",col.names=F,row.naes=F)
  
  load(file="data/results/betaInteraction/betaInteractionExIn.rda")
  
  variantAnnoNega <- read.delim("data/results/VEP/negative.delta.output.txt")
  variantAnnoPos <- read.delim("data/results/VEP/positive.delta.output.txt")
  variantAnnoNonSig <- read.delim("data/results/VEP/nonSig.delta.output.txt")
  
  rsID <- lapply(strsplit(as.character(res$rsID),";"),function(x){x})
  
#   rsIDENSG <- NULL
#   for (i in 1:nrow(eQTLPUTM))
#   {
#       
#     rsIDENSG <-  
#     
#     
#   }
  
  
  
  
  head(paste0(res$rsID,res$ge.gene))
  
  variantAnnoNonSig <- variantAnnoNonSig[which(paste0(variantAnnoNonSig$X.Uploaded_variation,variantAnnoNonSig$Gene) %in% paste0(res$rsID,res$ge.gene)),]
  ## remove duplicates
  variantAnnoNonSig <- unique(variantAnnoNonSig)
  
  ## get only the first consequence
  idx <- match(x = paste0(res$rsID,res$ge.gene), paste0(variantAnnoNonSig$X.Uploaded_variation,variantAnnoNonSig$Gene))
  
  idx<- idx[!is.na(idx)]
  
  variantAnnoNonSig <- variantAnnoNonSig[idx,]
  
  
  
  
  
  consNeg <- unlist(lapply(strsplit(as.character(variantAnnoNega$Consequence),","),function(x){x[1]}))
  consPos <- unlist(lapply(strsplit(as.character(variantAnnoPos$Consequence),","),function(x){x[1]}))
  consNonSig <- unlist(lapply(strsplit(as.character(variantAnnoNonSig$Consequence),","),function(x){x[1]}))
  
  rm(variantAnnoNega,variantAnnoNonSig,variantAnnoPos)
  nam <- names(sort(table(consNonSig),decreasing=T)/length(consNonSig))
    
  nonsig <- sort(table(consNonSig),decreasing=T)/length(consNonSig)
  positive <- sort(table(consPos),decreasing=T)/length(consPos)
  negative <- table(consNeg)/length(consNeg)
  
  counts <- rbind(nonsig=nonsig[nam],
                  positive=positive[nam],
                  negative=negative[nam])
  
  
  par(mar=c(15,5,3,2))
  barplot(counts, main="Variant Consequences",
          col=1:3,
          legend = rownames(counts), beside=TRUE,las=2,ylim=c(0,0.6))
  
  
  variantAnnoNega <- read.delim("data/results/VEP/negative.delta.output.txt")
  variantAnnoPos <- read.delim("data/results/VEP/positive.delta.output.txt")
  variantAnnoNonSig <- read.delim("data/results/VEP/nonSig.delta.output.txt")
  
  
  rm(variantAnnoNega,variantAnnoNonSig,variantAnnoPos)
  nam <- names(sort(table(consNonSig),decreasing=T)/length(consNonSig))
  
  nonsig <- sort(table(consNonSig),decreasing=T)
  positive <- sort(table(consPos),decreasing=T)
  negative <- table(consNeg)
  
  
  counts <- rbind(nonsig=nonsig[nam],
                  positive=positive[nam],
                  negative=negative[nam])
  
  
  counts[is.na(counts)]=0
  counts <- rbind(counts,total=apply(counts,2,sum))
  counts <- cbind(counts,total=apply(counts,1,sum))
  ftable(t(counts))
  ##counts <- t(counts) 
  
  
  ## we test the groups we have defined as non-sig,positive and negative
  pval <- NULL
  for(i in 1:(nrow(counts)-1)){
    
    pval[i] <- chisq.test(counts[c(4,i),1:ncol(counts)-1])$p.value
    
  }
    
  print(pval)
  rm(pval,i)
  
  ## we test the groups we have defined as non-sig,positive and negative but,
  ## instead of using the total number of sn for each categorie we use the sum of the two categories that were not tested
  pval <- NULL
  for(i in 1:(nrow(counts)-1)){
    
    pval[i] <- chisq.test(cbind(counts[4,1:ncol(counts)-1],
      counts[4,1:ncol(counts)-1]-counts[i,1:ncol(counts)-1]))$p.value
    
  }
  
  print(pval)
  rm(pval,i)
  
  
  
  ## with intronic as reference
  pval <- NULL
  for(i in 2:(ncol(counts)-1)){
    
    pval[i] <- chisq.test(counts[1:(nrow(counts)-1),c(1,i)])$p.value
    
  }
  
  print(pval)
  rm(pval,i)
  
  pval <- NULL
  for(i in 1:(ncol(counts)-1)){
    
    pval[i] <-chisq.test(counts[1:(nrow(counts)-1),c(15,i)])$p.value

  }
  
  print(pval)
  rm(pval,i)
  
  ## test for the non significant
  pval <- NULL
  for(i in 1:(ncol(counts)-1)){
    
    print(rbind(cbind(counts["nonsig",i],(counts[4,i]-counts["nonsig",i])),
                cbind(counts["nonsig",ncol(counts)],(counts[4,ncol(counts)]-counts["nonsig",ncol(counts)]))))
    pval[i] <-  chisq.test(rbind(cbind(counts["nonsig",i],(counts[4,i]-counts["nonsig",i])),
                                  cbind(counts["nonsig",ncol(counts)],(counts[4,ncol(counts)]-counts["nonsig",ncol(counts)]))))$p.value
    
    
  }
  
  print(pval)
  rm(pval,i)
  

  ## test for the positive
  pval <- NULL
  for(i in 1:(ncol(counts)-1)){
    
    print(rbind(cbind(counts["positive",i],(counts[4,i]-counts["positive",i])),
                cbind(counts["positive",ncol(counts)],(counts[4,ncol(counts)]-counts["positive",ncol(counts)]))))
    pval[i] <-  chisq.test(rbind(cbind(counts["positive",i],(counts[4,i]-counts["positive",i])),
                                 cbind(counts["positive",ncol(counts)],(counts[4,ncol(counts)]-counts["positive",ncol(counts)]))))$p.value
    
    
  }
  
  print(pval)
  rm(pval,i)
  
  
  
  ## test for the negative
  pval <- NULL
  for(i in 1:(ncol(counts)-1)){
    
    print(rbind(cbind(counts["negative",i],(counts[4,i]-counts["negative",i])),
                cbind(counts["negative",ncol(counts)],(counts[4,ncol(counts)]-counts["negative",ncol(counts)]))))
    pval[i] <-  chisq.test(rbind(cbind(counts["negative",i],(counts[4,i]-counts["negative",i])),
                                 cbind(counts["negative",ncol(counts)],(counts[4,ncol(counts)]-counts["negative",ncol(counts)]))))$p.value
    
    
  }
  
  print(pval)
  rm(pval,i)

}





## Examples of extremes for Mina
{ 

library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

load("data/results/betaInteraction/betaInteractionExIn.PUTM.rda")
tmp <-res 
load("data/results/betaInteraction/betaInteractionExIn.SNIG.rda")
res <- rbind(res,tmp)
rm(tmp)

res$FDRInter <- p.adjust(res$p.interaction,method="fdr",n=nrow(res))

## adding the one condition more for when the betas have oppositive signs

for(j in 1:nrow(res)){
  
  if(sign(res[j,"ge.beta"])!=sign(res[j,"gi.beta"])){
    res[j,"delta"] <- abs(res[j,"ge.beta"]) - abs(res[j,"gi.beta"])
  }else if(res[j,"ge.beta"] >= 0){
    res[j,"delta"] <- res[j,"ge.beta"] - res[j,"gi.beta"]
  }else
  {
    res[j,"delta"] <- res[j,"gi.beta"] - res[j,"ge.beta"]  
  }
}

par(mar=c(4, 4, 4, 4))
plot(res$delta,-log10(res$FDRInter),t="p",
     #     main=plot.title,
     xlab="delta",ylab="-log10(FDR)",main="Volcano plot FDR (PUTM+SNIG)",pch=21,cex=1)

##
fdr.threshold <- 0.05

for(j in 1:nrow(res)){
  
  if(res[j,"delta"] >= 0 & res[j,"FDRInter"] < fdr.threshold){
    res[j,"colors"] <- "red"
  }
  else if(res[j,"delta"] < 0 & res[j,"FDRInter"] < fdr.threshold){
    res[j,"colors"] <- "blue"  
  }else
    res[j,"colors"] <- "black"  
  res$p.interaction
}

plot(res$delta,-log10(res$FDRInter),t="p",col=res$colors,
     #     main=plot.title,
     xlab="Delta",ylab="-log10(FDR)",main="Volcano plot FDR (PUTM+SNIG)",bg=colors,pch=21,cex=1)
legend("topright",c("'positive'","'negative'","non significant"),col=c('red','blue','black'),pch=15)
rm(ensembl,j)

library(biomaRt)
## annotation of biotype
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","strand","gene_biotype"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=c(as.character(res$ge.gene)), mart=ensembl)



rm(ensembl)

rownames(geneNames) <- geneNames$ensembl_gene_id

res <- cbind(res,geneNames[as.character(res$ge.gene),c("external_gene_id","gene_biotype")])
rm(geneNames)

library(devtools)
library(R.utils)
setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

## POSTITIVE DELTAS
tmp <- res[which(res$delta > 0.6),] 

## NEGATIVE DELTAS
tmp <- res[which(res$delta< (-0.4)),] 

##
text(tmp$delta[5],-log10(tmp$FDRInter)[5],labels = tmp$external_gene_id[5],pos = 4)


gene <- "ENSG00000099290"
snp <- "chr10:51801386"

load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")
#exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
expr  <- resids[,as.character(gene)]
rm(IDs,resids)
rm(sampleInfo)



snp <- gsub("chr", "", snp)
load("data/general/imputed.dosage.rda")
load("data/general/imputed.info.rda")

dosage <- imputed.dosage[as.character(snp),names(expr)]
info <- imputed.info[which(as.character(imputed.info$marker) == as.character(snp)),]
table(round(as.numeric(dosage)))

identical(names(expr),names(dosage))
my.covTMP <- read.delim("data/general/eigenvec",sep=" ",header=F,row.names = 1)
head(my.covTMP)
my.cov0 <- as.matrix(my.covTMP[names(dosage),2:4])
eigenVectors <- t(my.cov0)
rm(my.cov0,my.covTMP)

title <- paste0("Gene Intronic ",tmp$external_gene_id[7]," (",snp,")")

par(mfrow=c(1,1))
plotEQTL(dosage =dosage,info = info,expr = expr,eigenVectors =eigenVectors,main = paste0("Gene exonic ",tmp$external_gene_id[7]," (",snp,")")  )


load("data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda")
#exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
expr  <- resids[,as.character(gene)]
rm(IDs,resids)
rm(sampleInfo)

plotEQTL(dosage =dosage,info = info,expr = expr,eigenVectors =eigenVectors,main = paste0("Gene intronic ",tmp$external_gene_id[7]," (",snp,")")  )


gene <- "ENSG00000224312"
snp <- "chr6:29833128"



load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")
#exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
expr  <- resids[,as.character(gene)]
rm(IDs,resids)
rm(sampleInfo)
snp <- gsub("chr", "", snp)

dosage <- imputed.dosage[as.character(snp),names(expr)]
info <- imputed.info[which(as.character(imputed.info$marker) == as.character(snp)),]
table(round(as.numeric(dosage)))

identical(names(expr),names(dosage))
my.covTMP <- read.delim("data/general/eigenvec",sep=" ",header=F,row.names = 1)
head(my.covTMP)
my.cov0 <- as.matrix(my.covTMP[names(dosage),2:4])
eigenVectors <- t(my.cov0)
rm(my.cov0,my.covTMP)


par(mfrow=c(1,2))
plotEQTL(dosage =dosage,info = info,expr = expr,eigenVectors =eigenVectors,main = paste0("Gene exonic ",tmp$external_gene_id[5]," (",snp,")")  )

load("data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda")
#exprFile <- "data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda"

load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
expr  <- resids[,as.character(gene)]
rm(IDs,resids)
rm(sampleInfo)

plotEQTL(dosage =dosage,info = info,expr = expr,eigenVectors =eigenVectors,main = paste0("Gene intronic ",tmp$external_gene_id[5]," (",snp,")")  )




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



## select significant and positive fold change greater than 2
idx <- which(res$FDRInter < 0.01 & res$delta >=0)             
positive <- res[idx,]  
#res <- res[-idx,]

## select significant and positive fold change lesser than -2
idx <- which(res$FDRInter < 0.01 & res$delta <0)             
negative <- res[idx,]  



par(mfrow=c(1,1))
hist(negative$DisGeneStart,col='skyblue',border=F,main= "TSS Intronic vs Exonic",
     sub=paste("Intronic:",length(negative$DisGeneStart),"Exonic:",length(positive$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(negative$DisGeneStart,positive$DisGeneStart)$p.value),freq=FALSE,breaks = 40)
hist(positive$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(negative$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(positive$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("Exonic","Intronic"),col=c('skyblue','red'),pch=15)

























