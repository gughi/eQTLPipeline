## with this script we get the closest gene to the intergenic region, basically we get the distances 

library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

getwd()

load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

## 5%FDR
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.05),1:7]
head(eQTLPUTM)

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
RPKM.cqn.PUTM <- RPKM.cqn
defReg.PUTM <- starStopReg
rm(PUTM,RPKM.cqn,covs)

eQTLPUTM <- cbind(eQTLPUTM,starStopReg[as.character(eQTLPUTM[,2]),])


library(GenomicRanges)

genesMap <- read.delim("data/general/ensemblGenes.txt")

head(genesMap)

## define the genomic regions
GR <- GRanges(seqnames = Rle(genesMap$Chromosome.Name),
              ranges = IRanges(start=genesMap$Gene.Start..bp.,
                               end = genesMap$Gene.End..bp.,
                               names = genesMap$Ensembl.Gene.ID))




plot(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre^2,las=2, main="correlation and distances (intergenic vs gene exonic quantification)",
     xlab="distance",ylab="r^2")
distanceNearest <- function(GRQuery,GR){
  
  x <- countOverlaps(GRQuery,GR,ignore.strand=T)
  idx <- nearest(GRQuery,GR,ignore.strand=T)
  dis <- distanceToNearest(GRQuery,GR,ignore.strand=T)
  nearGene <- names(GR[idx,])
  return(c(distance=as.data.frame(dis)$distance,gene=nearGene,overlap=x))
}



library(doParallel)
library(foreach)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(4)
clusterExport(cl, c("overlapsAny","countOverlaps","distanceToNearest","GRanges","IRanges","Rle","nearest"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
distance <- foreach(i=1:nrow(eQTLPUTM),.combine=rbind,.verbose=F)%dopar%distanceNearest(GRanges(seqnames = Rle(gsub("chr","",eQTLPUTM[i,"chr"])),
                                                                                           ranges = IRanges(start=eQTLPUTM[i,"start"],
                                                                                                            end = eQTLPUTM[i,"end"],
                                                                                                            names = genesMap[i,"gene"])),GR)
end <- Sys.time()
end-start
stopCluster(cl)
table(distance[,1]>"0")
rm(cl,GR,end,start)



hist(as.numeric(distance[which(as.numeric(distance[,1])>0),1]),breaks=40, main="histogram of distance btw intergenic and nearest gene", 
     xlab=" distance in bp")

hist(as.numeric(distance[which(as.numeric(distance[,1])>0 & as.numeric(distance[,1])<10000),1]),breaks=40, main="histogram of distance btw intergenic and nearest gene", 
     xlab=" distance in bp")

table(as.numeric(distance[,1])>0)
# FALSE  TRUE 
# 878  1256 


table(distance[,3])


# load("data/expr/rawCounts/intergenic/transcribedRegions/PUTM/allChromosomes.rda")
# ## load data has 2 different objects: annotation,coverage
# ## load the sample info
# load("data/general/sampleInfo.rda")
# ## we select only samples from PUTM
# 
# PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
# rm(sampleInfo)
# 
# exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
# colnames(exprIntergenic) <- c("chr","start","end","width",colnames(coverage))
# rm(annotation,coverage)    
# 
# 
# ## select regions with length => 100bp
# exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
# # We don't do filtering since it doesn't make any since for DERFINDER; 
# # regions are detected based on the expression
# # creation of identifiers for the regions
# IDs <- paste0("DER",c(1:nrow(exprIntergenic)))
# rownames(exprIntergenic) <- IDs
# rm(IDs)
# 
# save(exprIntergenic,file="data/expr/rawCounts/intergenic/exprIntergenic.PUTM.rda") 



eQTLPUTM <- cbind(eQTLPUTM,distance)
colnames(eQTLPUTM)[2] <- "regionID"

intExpr <- RPKM.cqn.PUTM
rm(RPKM.cqn.PUTM)
load("data/expr/rawCounts/genic/exprSQ.rda")

load("data/expr/rawCounts/intergenic/exprIntergenic.PUTM.rda")


load("data/general/sampleInfo.rda")
## we select only samples from PUTM
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]  

correla <- NULL
for (i in 1:nrow(eQTLPUTM))
{
  correla[i] <-  cor(as.numeric(exprIntergenic[as.character(eQTLPUTM$regionID[i]),as.character(PUTM$A.CEL_file)]),
              as.numeric(exprSQ[as.character(eQTLPUTM$gene[i]),as.character(PUTM$A.CEL_file)]))
  
}


eQTLPUTM <- cbind(eQTLPUTM,corre=correla)

eQTLPUTMTmp <- eQTLPUTM[order(eQTLPUTM$distance),]

## removing regions that have correlation equal to 0 to visualisation perposes 
eQTLPUTMTmp <- eQTLPUTMTmp[-which(is.na(eQTLPUTMTmp$corre)),]

plot(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre^2,las=2, main="correlation and distances (intergenic vs gene exonic quantification)",
     xlab="distance",ylab="r^2")

cor(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre)

##eQTLPUTMTmp <- eQTLPUTMTmp[-which(as.numeric(as.character(eQTLPUTMTmp$distance))==0),]

plot(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre^2,las=2, main="correlation and distances (intergenic vs gene exonic quantification)",
     xlab="distance",ylab="r^2")


abline(glm((eQTLPUTMTmp$corre^2)~as.numeric(as.character(eQTLPUTMTmp$distance))),col="red")

abline(v=10000,col="blue")
abline(h=0.1,col="green")

plot(as.numeric(as.character(eQTLPUTMTmp$width)),eQTLPUTMTmp$corre^2,las=2, main="width and correlation (intergenic vs gene exonic quantification)",
     xlab="width",ylab="r^2")

abline(glm((eQTLPUTMTmp$corre^2)~as.numeric(as.character(eQTLPUTMTmp$width))),col="red")




load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")

head(eQTLPUTMTmp)


plot(as.numeric(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"distance"]),eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"]
     ,las=2, main="correlation and distances (intergenic vs gene exonic quantification)",
     xlab="distance",ylab="r^2")


plot(as.numeric(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"width"]),eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"]
     ,las=2, main="correlation and distances (intergenic vs gene exonic quantification)",
     xlab="distance",ylab="r^2")


barplot(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"])

t.test(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"],
       eQTLPUTMTmp[-which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"])

mean(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"])
mean(eQTLPUTMTmp[-which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"])


# plot(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre^2,las=2, main="correlation and distances (intergenic vs gene exonic quantification)",
#      xlab="distance",ylab="r^2")
# points(eQTLPUTMTmp[-which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"distance"],(eQTLPUTMTmp[-which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"]^2),col="red")
# points(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"distance"],(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),"corre"]^2),col="blue")


barplot(as.numeric(eQTLPUTMTmp[which(eQTLPUTMTmp$gene %in% eQTLPUTM$gene),]))


length(unique(as.character(eQTLPUTMTmp$gene)))










