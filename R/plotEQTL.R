plotEQTL <- function(dosage,info,expr,eigenVectors,main=NULL)
{
  ##snp <- "chr6:29955809"
  # dosageFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.dosage"
  # infoFile <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"
  
## load the expresssion
  
  stopifnot(identical(names(expr),colnames(eigenVectors)))
  stopifnot(identical(names(expr),names(dosage)))
  
  
  plot(round(as.numeric(dosage)),expr,xaxt="n",ylab="expression",xlab="Genotyope", main=main)
  abline(glm(expr ~ as.numeric(dosage)+ eigenVectors[1,]+eigenVectors[2,]+eigenVectors[3,]),col="red")
  fit <- coef(summary(glm( expr ~ as.numeric(dosage) + eigenVectors[1,]+eigenVectors[2,]+eigenVectors[3,]) ))
  
  mtext(paste("pval",fit["as.numeric(dosage)", "Pr(>|t|)"]) , side=1, at=1, cex=0.6, font=2, line=2 )  ## add tissue and pvals
  
  mtext( side=1, paste0(as.character(info[1,"Al2"]), as.character(info[1,"Al2"])), at=(0), cex=0.6, line=-0.25 )
  mtext( side=1, paste0(as.character(info[1,"Al2"]), as.character(info[1,"Al1"])), at=(1),   cex=0.6, line=-0.25)
  mtext( side=1, paste0(as.character(info[1,"Al1"]), as.character(info[1,"Al1"])), at=(2), cex=0.6, line=-0.25 )
}



