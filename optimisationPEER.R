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


##for (i in 1:20) 
##{
##  doPEER(RPKM.cqn=t(RPKM.cqn),nFactors=13,
##         covs=covs[as.character(colnames(RPKM.cqn)),c("Age","Gender","Region")]
##         ,outputFile=paste0("testPEER/RNDMPEERrun",i))
##}


## this is actually the library I executed 
 library(doParallel)
 library(foreach)
  cl <- makeCluster(20)
  registerDoParallel(cl)
  foreach(i=1:20, .verbose=T)%dopar%system(paste0("R --vanilla --file=/home/seb/projectsR/eQTLPipeline/optimisationPEER.R --args ",i," > /home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i,".log"))
  stopCluster(cl)
  rm(cl)
## test
## foreach(i=1:20, .verbose=T)%dopar%system(paste0("echo ",i," > /home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i,".log"))


# correPlot <- function (mat1,mat2,xlab)
# {
#  library(gplots)
#   linp<-matrix(ncol=ncol(mat1),nrow=ncol(mat2))  
#   rownames(linp)<-colnames(mat2)
#   colnames(linp)<-colnames(mat1)
#   rsquared<-matrix(ncol=ncol(mat1),nrow=ncol(mat2))
#   rownames(rsquared)<-colnames(mat2)
#   colnames(rsquared)<-colnames(mat1)
#   for (i in 1:ncol(mat2)){
#     for (j in 1:ncol(mat1)){
#       fit<-lm(mat1[,j]~mat2[,i])
#       s<-summary(fit)
#       linp[i,j]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
#       rsquared[i,j]<-s$r.squared[1]
#     }}
#   
#   
#   smallest=-20
#   linp10<-log10(linp)
#   linp10<-replace(linp10,linp10<=smallest,smallest)
#   print(linp10)
#   
#   rsquaredTMP <- rsquared
#   rsquaredTMP[which(linp10>-5)] = NA 
#   heatmap.2(linp10,Colv=F,Rowv=F,dendrogram="none",trace="none",symbreaks=F,symkey=F,breaks=seq(-20,0,length.out=100),key=T,col=heat.colors(99),
#             cexRow=1,cexCol=1,colsep=NULL,rowsep=NULL,sepcolor=sepcolor,sepwidth=sepwidth,
#             main="",labCol=paste(1:ncol(linp10),sep=""),margins=c(5,7),labRow=,xlab=xlab,
#             cellnote=rsquaredTMP,notecol="black",notecex=1) 
#   
# }
# 
# 
# load("data/general/RPKMCQNcovs.rda")
# 
# corPEER <- data.frame(matrix(NA, nrow=15, ncol=20))
# for (i in 1:20) 
#   {
#           PEER <- read.csv(paste0("/home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i), row.names=1)
#           # do plot
#            jpeg(paste0("/home/seb/projectsR/eQTLPipeline/testPEER/RNDMPEER",i,".jpeg"))
#            correPlot(PEER,covs[rownames(PEER),],paste("PEER random",i,"vs Known factors"))
#            dev.off()
#           # calculate correaltion
#           corPEER[,i] <- as.data.frame(apply(PEER[,c(1:2,4:16)],2,function(x) cor(PEER[,3],x)))      
#           colnames(corPEER)[i] <- paste0("RNDMPEER",i)
#           rm(PEER)
#   }
#   rm(i)
# 
# save(corPEER,file="testPEER/correlation.rda")
# 
# rm(RPKM.cqn)
# rm(cortmp)
# rm(covs)
# 
