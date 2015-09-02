#' Function that does that returns the residuals 
#' @param exprThe expression matrix nxm where n are the samples and m the genes
#' @param covariates to use for the correction
#' @param covs the covariates to exclude from PEER analysis
#' @param outputFile where the results are going to be saved 
doResidualCorrection <-
function(expr,covs,outputFile)
{
  
  ## load the knowing factors (age and sex ) and unknown factors (PEER axes)
  ## Residual correction with linear model, step 3
  print("Performing the residual correction")
  stopifnot(identical(rownames(expr),rownames(covs)) )
  resids <- apply(expr, 2, function(y){
    lm( y ~ . , data=covs)$residuals 
  })
  
  print("Saving the residuals")
  save(resids,file=paste0(outputFile))
  resids
}
