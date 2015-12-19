# beta intercation test between intronic and exonic eQTLs


## calculate the beta interaction for my eQTLs
## the FUNCTION TO CALCULATE THE INTERACTION IS AT THE BOTTOM
## 27-1-14


library(lme4)



## load the eQTLs
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLExonic <- eQTLPUTM[,c("snps","gene")]
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLIntronic <- eQTLPUTM[,c("snps","gene")]
rm(eQTLPUTM)

## get the beta interaction for intronic and exonic
eQTLExonic$type <- "exonic"
eQTLIntronic$type <- "intronic"

eQTLs <- rbind(eQTLExonic,eQTLIntronic)

## load the expression(residual corrected values)
load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")
exprExonic <- resids
load("data/expr/normalisedCounts/genic/geneIntronic/resids.PUTM.rda")
exprIntronic <- resids
rm(resids)


for (j in 1:nrow(eQTLs))
{
  ## load(paste0("/home/seb/forMar/RemovedBatchAgeGenderRIN/expr/cases/t",eQTLcases$tID[j],".rda"))
  ## eQTLs$gene
  j<-1
  
  ## load the expression for each individual gene for exonic and intronic
  exprE <- exprExonic[,as.character(eQTLs$gene[j])]
  exprI <- exprIntronic[,as.character(eQTLs$gene[j])]
  
  ## add to the expression the tissue information or the exonic/intronic information
  expr <- cbind(exprE,"Exonic")
  exprtmp <- cbind(exprI,"Intronic")
  expr <- rbind(expr,exprtmp)
  rm(exprtmp,exprE,exprI)
  
  ## merge the ctrl samples
  #load(paste0("/home/seb/forMar/RemovedBatchAgeGenderRIN/expr/ctrl/t",eQTLcases$tID[j],".rda"))
  #tExprCases <- t(tExprCases)
  #expr0tmp <- data.frame(tExprCases[as.character(eQTLcases$exprID[j]),])
  #colnames(expr0tmp) <- as.character(eQTLcases$exprID[j])
  #expr <- expr0tmp ## to calculate the p value of the controls
  #expr0tmp$condition <- "ctrl"
  #expr0 <- rbind(expr0,expr0tmp)
  #rm(expr0tmp)
  ## expr0 <- expr0tmp
  
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
  
  head(sampleInfo)
  head(rownames(expr))
  
  ##if (eQTLs$tissue[j] == "PUTM")
  ##{
  ##  selSamples <- sampleInfo[which(sampleInfo$A.CEL_file %in% rownames(exprSNIG)),1]
  ##  snps <- as.numeric(markers[as.character(eQTLs$snps[j]),selSamples])
  ##  snps <- round(snps)
  ##  exprtmp <- exprSNIG[,eQTLs$gene[j]]
  ##  store <- t( apply( as.matrix(exprtmp), 2, function(e){
  ##    fit <- coef( summary( glm( e ~ snps ) ) )
  ##    return( c( fit[ , "Estimate"], fit["snps", "Pr(>|t|)"] ) )
  ##  }) )
  ##  rm(snps)
  ##}
  ##else
  ##{
  ##  selSamples <- sampleInfo[which(sampleInfo$A.CEL_file %in% rownames(exprPUTM)),1]
  ##  snps <- as.numeric(markers[as.character(eQTLs$snps[j]),selSamples])
  ##  snps <- round(snps)
  ##  exprtmp <- exprPUTM[,eQTLs$gene[j]]
  ##  store <- t( apply( as.matrix(exprtmp), 2, function(e){
  ##    fit <- coef( summary( glm( e ~ snps ) ) )
  ##    return( c( fit[ , "Estimate"], fit["snps", "Pr(>|t|)"] ) )
  ##  }) )
  ##  rm(snps)
  ##}
  colnames(snp0) <- rownames(expr)
    
  f <- Beta.interaction(expr,snp0)
  restmp <- cbind(store,f)
  colnames(restmp) <- c("Intercept.Beta", "Beta.ctrl", "Pvalue.ctrl","p.interaction")
  
  if (nrow(res)==0)
  {
    res <- restmp
  } else {
    res <- rbind(res,restmp)
  }
  rm(store,f)
  
  
}

completeQTLs <- cbind(eQTLcases,res)

completeQTLs$Intercept.Beta <- NULL
write.csv(completeQTLs,file="/home/seb/forMar/eQTLs/RemovedBatchAgeGenderRIN/sentinalised1FDR/completeEQTL.cases.csv")

eQTLcases <- read.table("/home/seb/forMar/eQTLs/RemovedBatchAgeGenderRIN/sentinalised1FDR/mTLE1FDRAnnotated.cases.tsv")
res<-data.frame()

## for controls

eQTLctrl <- read.table("/home/seb/forMar/eQTLs/RemovedBatchAgeGenderRIN/sentinalised1FDR/mTLE1FDRAnnotated.ctrl.tsv")
res<-data.frame()

for (j in 1:nrow(eQTLctrl))
{
  load(paste0("/home/seb/forMar/RemovedBatchAgeGenderRIN/expr/ctrl/t",eQTLctrl$tID[j],".rda"))
  tExprCases <- t(tExprCases)
  expr0 <- data.frame(tExprCases[as.character(eQTLctrl$exprID[j]),])
  colnames(expr0) <- as.character(eQTLctrl$exprID[j])
  expr0$condition <- "ctrl"
  
  ## merge the ctrl samples
  load(paste0("/home/seb/forMar/RemovedBatchAgeGenderRIN/expr/cases/t",eQTLctrl$tID[j],".rda"))
  tExprCases <- t(tExprCases)
  expr0tmp <- data.frame(tExprCases[as.character(eQTLctrl$exprID[j]),])
  colnames(expr0tmp) <- as.character(eQTLctrl$exprID[j])
  expr <- expr0tmp ## to calculate the p value of the controls
  expr0tmp$condition <- "cases"
  expr0 <- rbind(expr0,expr0tmp)
  rm(expr0tmp)
  ## expr0 <- expr0tmp
  
  
  load(paste0("/home/seb/forMar/snpsInfo/byTranscript/t",eQTLctrl$tID[j],".rda"))
  snp0 <- snps[as.character(eQTLctrl$snps[j]),rownames(expr0)]
  snp0 <- round(snp0)
  
  
  snps <- as.numeric(snps[as.character(eQTLctrl$snps[j]),rownames(expr)])
  snps <- round(snps)
  store <- t( apply( expr, 2, function(e){
    fit <- coef( summary( glm( e ~ snps ) ) )
    return( c( fit[ , "Estimate"], fit["snps", "Pr(>|t|)"] ) )
  }) )
  rm(snps)
  
  f <- Beta.interaction(expr0,snp0)
  restmp <- cbind(store,f)
  colnames(restmp) <- c("Intercept.Beta", "Beta.cases", "Pvalue.cases","p.interaction")
  
  if (nrow(res)==0)
  {
    res <- restmp
  } else {
    res <- rbind(res,restmp)
  }
  rm(store,f)
  
  
}

completeQTLs <- cbind(eQTLctrl,res)

completeQTLs$Intercept.Beta <- NULL
write.csv(completeQTLs,file="/home/seb/forMar/eQTLs/RemovedBatchAgeGenderRIN/sentinalised1FDR/completeEQTL.ctrl.csv")





Beta.interaction <- function(expr0, snp0){
  
  
  ## Check whether the expr and the snps have the same order of samples
  stopifnot( identical( rownames(expr0), colnames(snp0) ) )
  ## stop if number of snp is different to 1
  stopifnot( nrow(snp0) == 1 )
  
  ## replicate the values of snp0 for the number of probes in expr0
  my.snp  <- rep( as.numeric(snp0), (ncol(expr0)-1) )
  
  ## make the expression in a single vector
  my.expr <- c(t(as.matrix(expr0[,1]))) # vectorize it
  
  ## give a index for each tissue of exon
  ## my.exon <- factor( rep( 1:nrow(expr0), each=ncol(expr0) ) )  ## exon/tissue number
  ## my exon is actually the condition
  my.exon <- factor( expr0[,2] )
  
  ## give index to the individuals
  my.ind  <- factor( rep( 1:nrow(expr0), ncol(expr0) ) ) ## individuals number
  
  ## create data frame aggregating per colonne
  df <- cbind.data.frame(my.expr, my.exon, my.snp, my.ind)
  ## select by complete cases, look for complete.cases
  df <- df[ complete.cases(df), ]
  ## convert as factors column my.exon
  df$my.exon <- as.factor(df$my.exon)
  ## convert as factors column my.ind
  df$my.ind <- factor(df$my.ind)
  ## remove all the object a part of df
  rm(my.expr, my.exon, my.snp, my.ind, expr0, snp0)
  
  ## mixed model approach, fixed sQTL effects (MW recommended)
  mod0 <- lmer( my.expr ~  1 + (1|my.ind) + (1|my.exon) + my.snp , data=df )
  mod1 <- lmer( my.expr ~  1 + (1|my.ind) + (1|my.exon) + my.snp + my.exon:my.snp , data=df )
  p <- anova(mod0, mod1)$"Pr(>Chisq)"[2]
  rm(mod0, mod1, df)
  return( p )
}




