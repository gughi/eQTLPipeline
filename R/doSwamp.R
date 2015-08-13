doSwamp <-
function(RPKM.cqn,covs)
{
  library(swamp)
  RPKM.cqn.tmp <- t(RPKM.cqn)
  expr.data.o <- RPKM.cqn.tmp[order(rownames(RPKM.cqn.tmp)),]
  
  ## Order your data in the same way
  traits.o <- covs[as.character(rownames(expr.data.o)),]
  rownames(traits.o) <- rownames(expr.data.o)
  
  
  print(head(covs))
  stopifnot(identical(rownames(traits.o),rownames(expr.data.o)))
  
  
  #Do PCA analysis with up to 10 axis
  res1 <- prince(t(expr.data.o),traits.o,top=15,permute=TRUE)
  #PLOT WHAT YOU NEED
  prince.plot(prince=res1)
  
  #Plot the variation included by each PCA axis
  res2 <- prince.var.plot(t(expr.data.o),show.top=50,npermute=10)
  
  hca.plot(as.matrix(t(expr.data.o)),as.data.frame(traits.o))

}
