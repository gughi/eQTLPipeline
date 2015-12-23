## beta interaction test
betaInteraction <- function(expr0, snp0){
  
  
  ## Check whether the expr and the snps have the same order of samples
  stopifnot(identical( rownames(expr0), colnames(snp0)))
  ## stop if number of snp is different to 1
  stopifnot( nrow(snp0) == 1)
  
  ## replicate the values of snp0 for the number of probes in expr0
  my.snp  <- rep( as.numeric(snp0),(ncol(expr0)-1))
  
  ## make the expression in a single vector
  my.expr <- c(t(as.matrix(expr0[,1]))) # vectorize it
  
  ## give a index for each tissue of exon
  ## my.exon <- factor( rep( 1:nrow(expr0), each=ncol(expr0) ) )  ## exon/tissue number
  ## my exon is actually the condition
  my.exon <- factor(expr0[,2])
  
  ## give index to the individuals
  ##my.ind  <- factor(rep( 1:ncol(expr0),nrow(expr0)-1)) ## individuals number
  my.ind <- as.numeric(as.factor(rownames(expr0)))
  
  
  ## create data frame aggregating per colonne
  df <- cbind.data.frame(my.expr, my.exon, my.snp, my.ind)
  ## select by complete cases, look for complete.cases
  df <- df[ complete.cases(df), ]
  ## convert as factors column my.exon
  df$my.exon <- as.numeric(factor(df$my.exon))
  ## convert as factors column my.ind
  df$my.ind <- as.numeric(factor(df$my.ind))
  ## remove all the object a part of df
  df$my.expr <- as.numeric(as.character(df$my.expr))
  df$my.snp <- as.numeric(as.character(df$my.snp))
  rm(my.expr, my.exon, my.snp, my.ind, expr0, snp0)
  
  
  ## mixed model approach, fixed sQTL effects (MW recommended)
  mod0 <- lmer( my.expr ~  1 + (1|my.ind) + my.exon + my.snp , data=df )
  mod1 <- lmer( my.expr ~  1 + (1|my.ind) + my.exon + my.snp + my.exon:my.snp , data=df )
  p <- anova(mod0, mod1)$"Pr(>Chisq)"[2]
  rm(mod0, mod1, df)
  return( p )
}



## beta interaction test
betaInteractionExons <- function(expr0, snp0){
    
  ## Check whether the expr and the snps have the same order of samples
  stopifnot(identical( rownames(expr0), colnames(snp0)))
  ## stop if number of snp is different to 1
  stopifnot( nrow(snp0) == 1)
  
  ## replicate the values of snp0 for the number of probes in expr0
  my.snp  <- rep( as.numeric(snp0),(ncol(expr0)-1))
  
  ## make the expression in a single vector
  my.expr <- c(t(as.matrix(expr0[,1]))) # vectorize it
  
  ## give a index for each tissue of exon
  ## my.exon <- factor( rep( 1:nrow(expr0), each=ncol(expr0) ) )  ## exon/tissue number
  ## my exon is actually the condition
  my.exon <- factor(expr0[,2])
  
  ## give index to the individuals
  ##my.ind  <- factor(rep( 1:ncol(expr0),nrow(expr0)-1)) ## individuals number
  my.ind <- as.numeric(as.factor(rownames(expr0)))
  
  
  ## create data frame aggregating per colonne
  df <- cbind.data.frame(my.expr, my.exon, my.snp, my.ind)
  ## select by complete cases, look for complete.cases
  df <- df[ complete.cases(df), ]
  ## convert as factors column my.exon
  df$my.exon <- as.numeric(factor(df$my.exon))
  ## convert as factors column my.ind
  df$my.ind <- as.numeric(factor(df$my.ind))
  ## remove all the object a part of df
  df$my.expr <- as.numeric(as.character(df$my.expr))
  df$my.snp <- as.numeric(as.character(df$my.snp))
  rm(my.expr, my.exon, my.snp, my.ind, expr0, snp0)
  

  ## mixed model approach, fixed sQTL effects (MW recommended)
  mod0 <- lmer( my.expr ~  1 + (1|my.ind) + (1|my.exon) + my.snp , data=df )
  mod1 <- lmer( my.expr ~  1 + (1|my.ind) + (1|my.exon) + my.snp + my.exon:my.snp , data=df )
  p <- anova(mod0, mod1)$"Pr(>Chisq)"[2]
  rm(mod0, mod1, df)
  return( p )
}



betaInteractionExonicIntronic <- function(expr0, snp0){
  
  
  ## Check whether the expr and the snps have the same order of samples
  stopifnot(identical( rownames(expr0), colnames(snp0)))
  ## stop if number of snp is different to 1
  stopifnot( nrow(snp0) == 1)
  
  ## replicate the values of snp0 for the number of probes in expr0
  my.snp  <- rep( as.numeric(snp0),(ncol(expr0)-1))
  
  ## make the expression in a single vector
  my.expr <- c(t(as.matrix(expr0[,1]))) # vectorize it
  
  ## give a index for each tissue of exon
  ## my.exon <- factor( rep( 1:nrow(expr0), each=ncol(expr0) ) )  ## exon/tissue number
  ## my exon is actually the condition
  my.exon <- factor(expr0[,2])
  
  ## give index to the individuals
  ##my.ind  <- factor(rep( 1:ncol(expr0),nrow(expr0)-1)) ## individuals number
  my.ind <- as.numeric(as.factor(rownames(expr0)))
  
  ## get the genetic PCAs
  genPCA <- getGeneticPCA(unique(rownames(expr0)))
  genPCA <- genPCA[rownames(expr0),]
  
  
  ## create data frame aggregating per colonne
  df <- cbind.data.frame(my.expr, my.exon, my.snp, my.ind)
  ## select by complete cases, look for complete.cases
  df <- df[ complete.cases(df), ]
  ## convert as factors column my.exon
  
  df$my.exon <- as.numeric(factor(df$my.exon))
  df$my.intronic <- df$my.exon -1
  df$my.exonic <- 1-df$my.intronic
  ## convert as factors column my.ind
  df$my.ind <- as.numeric(factor(df$my.ind))
  ## remove all the object a part of df
  df$my.expr <- as.numeric(as.character(df$my.expr))
  df$my.snp <- as.numeric(as.character(df$my.snp))
  rm(my.expr, my.exon, my.snp, my.ind, expr0, snp0)
  
  
  fit    <- glm( my.expr ~ my.snp*my.exonic + genPCA[,2]+genPCA[,3]+genPCA[,4],data=df)
  fitRev <- glm( my.expr ~ my.snp*my.intronic+ genPCA[,2]+genPCA[,3]+genPCA[,4],data=df )
  keepCols <- c("Estimate", "Std. Error", "Pr(>|t|)")
  p <- c(coef(summary(fit))[ "my.snp:my.exonic", keepCols ])
  
  rm(fit, fitRev, df)
  return( p )
}


getGeneticPCA <- function(sampleIDs)
{
  eigenFile <- "/home/seb/plinkOutput/eigenvec"
  ## load the samples information
  load("data/general/sampleInfo.rda")  
  IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(sampleIDs)),"U.SD_No"]
  IDs <- gsub("/","_",IDs)
  
  my.covTMP <- read.table.rows(paste0(eigenFile), keepRows=IDs, sep=" ",header=F)
  my.covTMP <- my.covTMP[as.character(IDs),]
  rownames(my.covTMP) <- sampleIDs
  return(my.covTMP)
}

