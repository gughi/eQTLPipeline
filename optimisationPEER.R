## optimisation PEER


#' Function that calculates hidden factors in the experiment usin PEER 
#' @param RPKM.cqn = The expression normalised
#' @param nFactors the number of hidden factor we want to discover
#' @param covs the covariates to exclude from PEER analysis
#' @param outputFile where the results are going to be saved
print( i <- as.numeric( commandArgs(trailingOnly=T)[1] ) )

doPEER <- function(RPKM.cqn,nFactors,covs,outputFile)
{
  
  library(peer,lib="/home/seb/RlibraryTest/")
  
  
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
  
  PEER_setVarTolerance(model, 0.0000000000000000000000001)
  
  
  ## run PEER, step 5
  print(paste("variance tolerance",PEER_getVarTolerance(model)))
  PEER_update(model)
  
  ## get the hidden factors
  factor <- PEER_getX(model)
  head(factor)
  rownames(factor) <- rownames(RPKM.cqn)
  
  precision <- PEER_getAlpha(model)
  print(dim(precision))
  ##plot(precision)
  
  write.csv(factor,file=outputFile)
  
  
}

load("/home/seb/projectsR/eQTLPipeline/data/general/RPKMCQNcovs.rda")

doPEER(RPKM.cqn=t(RPKM.cqn),nFactors=13,
       covs=covs[as.character(colnames(RPKM.cqn)),c("Age","Gender","Region")]
       ,outputFile=paste0("/home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i))

print(" PEER Executed ")


# # this is actually the library I executed 
#  library(doParallel)
#  library(foreach)
#   cl <- makeCluster(20)
#   registerDoParallel(cl)
#   foreach(i=1:20, .verbose=T)%dopar%system(paste0("R --vanilla --file=/home/seb/projectsR/eQTLPipeline/optimisationPEER.R --args ",i," > /home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i,".log"))
#   stopCluster(cl)
#   rm(cl)
## test
## foreach(i=1:20, .verbose=T)%dopar%system(paste0("echo ",i," > /home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i,".log"))
# 
# 
# 
# 
#  load("/home/seb/projectsR/eQTLPipeline/data/general/RPKMCQNcovs.rda")
#  
#  corPEER <- data.frame(matrix(NA, nrow=15, ncol=20))
#  for (i in 1:20) 
#    {
#            PEER <- read.csv(paste0("/home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i), row.names=1)
#            # do plot
#             jpeg(paste0("/home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i,".jpeg"))
#             correPlot(PEER,covs[rownames(PEER),],paste("PEER random",i,"vs Known factors"))
#            dev.off()
#           # calculate correaltion
#           corPEER[,i] <- as.data.frame(apply(PEER[,c(1:2,4:16)],2,function(x) cor(PEER[,3],x)))      
#           colnames(corPEER)[i] <- paste0("RNDMPEER",i)
#           rm(PEER)
#   }
#   rm(i)
# 
# save(corPEER,file="/home/seb/projectsR/eQTLPipeline/testPEER/correlation.rda")
# 
# rm(RPKM.cqn)
# rm(cortmp)
# rm(covs)
# 

# How I select the best PEER I basically calculate the rsquared for all the PEER test,
# from each test calculate the maximum and then I select the PEER test that has minimum rsquared beetween all the tests  
# load("testPEER/correlation.rda")
# print(paste("Test with minimun correlation",match(min(apply((corPEER^2),2,max)),apply((corPEER^2),2,max))))
# 
# load("data/general/RPKMCQNcovs.rda")
# PEERRNDPEER20 <- read.csv("testPEER/RNDMPEER20",row.names=1)
# correPlot(PEERRNDPEER20,covs[rownames(PEERRNDPEER20),],"Correlation 'best' PEER and known factors")
# 
# 
# cor(covs$Region,covs$uniqueMappedRead)
#
#

