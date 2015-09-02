#' Function that calculates hidden factors in the experiment usin PEER 
#' @param RPKM.cqn = The expression normalised
#' @param nFactors the number of hidden factor we want to discover
#' @param covs the covariates to exclude from PEER analysis
#' @param outputFile where the results are going to be saved
doPEER <-
function(RPKM.cqn,nFactors,covs,outputFile)
{
  
  library(peer,lib="/home/seb/Rlibrary/")
  
  
  ## WE DON'T INCLUDE THE GENETIC COVARIATS
  ## geneticPCA <- read.table.rows("/home/guelfi/plinkOutput/eigenvec", keepRows=covs$U.SD_No, sep=" ",header=F)
  ## geneticPCA <- geneticPCA[covs$U.SD_No,]
  ## rownames(geneticPCA) <- rownames(covs)
  
  ## WE INCLUDE AGE,GENDER AND TISSUE
  covs <- covs[as.vector(rownames(RPKM.cqn)),]
  stopifnot(identical(rownames(covs),rownames(RPKM.cqn)))
  
  ## create the PEER model
  model=PEER()
  ## set to obtain only the first 15 hidden factors
  PEER_setNk(model, nFactors)
  ## set the covariants, step 4
  PEER_setCovariates(model, data.matrix(covs))
  PEER_setPhenoMean(model,as.matrix(RPKM.cqn))
  
  
  ## check if the covariants have been updatef properly
  print(paste("dimension of the expression matrix",print(dim(PEER_getPhenoMean(model)))))
  print(paste("dimension of the covariate matrix",print(dim(PEER_getCovariates(model)))))
  ##head(PEER_getPhenoMean(model))
  
  ## run PEER, step 5
  PEER_update(model)
  
  ## get the hidden factors
  factor <- PEER_getX(model)
  head(factor)
  rownames(factor) <- rownames(RPKM.cqn)
  
  write.csv(factor,file=outputFile)
}
