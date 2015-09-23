# Differentially expressed analysis

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade") 
# biocLite("DESeq2")
# biocLite("easyRNASeq")
# biocLite("cqn")


library(DESeq2)

setwd("/home/guelfi/eQTLPipeline/")
load("data/expr/rawCounts/genic/exprSQ.rda")

load("data/general/sampleInfo.rda")

exprSQ[is.na(exprSQ)]=0
## remove genes that not expressed in any gene
exprSQ <- exprSQ[rowSums(exprSQ>0)>0,]

countData <- exprSQ[,as.character(sampleInfo$A.CEL_file)] 
rm(exprSQ)
colData <- sampleInfo[,c("U.Region_simplified","A.Ovation.Prep.Batch.of.succesfull.prep")]
colnames(colData) <- c("region","batch")
rownames(colData) <- sampleInfo$A.CEL_file

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData = colData,
                              design =  ~  region)


#####################################################################
## now we calculate the CQN offset to correct for GC content bias ###
#####################################################################

library(devtools)
load_all()

library(doParallel)
library(foreach)

exonsdef <- read.csv("data/general/exonDef.csv")


## calculation of genes only exons length

library(doParallel)
library(foreach)

nCores <- 7
detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(nCores)
clusterExport(cl, c("getRegionsWidth","defExonicRegions"))

registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
geneswidth <- foreach(i=1:length(rownames(exprSQ)),.export=c("getRegionsWidth","defExonicRegions"),.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(exprSQ)[i],exonsdef)
#geneswidth <- foreach(i=1:10,.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(expr)[i],exonsdef)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)

length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[as.character(rownames(exprSQ))]


cl <- makeCluster(nCores)
clusterExport(cl, c("getRegionsBED","defExonicRegions"))

registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
cat(paste("calculating GC content...","\n"))
exonicRegions <- foreach(i=1:length(rownames(exprSQ)),.combine=rbind,.verbose=F)%dopar%getRegionsBED(rownames(exprSQ)[i],exonsdef)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl)

write.table(data.frame(exonicRegions), file = paste0("data/general/exonicRegions.BED"), row.names = F, 
            col.names = F, quote = F)

## we filter things that not match with the fasta file

system("grep -v HG* data/general/exonicRegions.BED  | grep -v LRG* | grep -v HS* | cat > data/general/exonicRegionsFiltered.BED ")

cmd <- paste0("/apps/BEDTools/2.24.0/bin/bedtools nuc -fi /home/ukbec/bowtie2Index/genome37.72.fa -bed data/general/exonicRegionsFiltered.BED > data/general/GCcontRegionsExonic")

## calculate GC content with bedtools
system(cmd)

end <- Sys.time()
end-start
cat(paste("GC content calculated in",end-start,"\n"))
rm(end,start)
GCcontentTab <- read.delim("data/general/GCcontRegionsExonic")
cat(paste("GC content saved in data/general/GCcontRegionsExonic","\n"))

rm(cmd)
## detectCores()
## [1] 24

cl <- makeCluster(nCores)
clusterExport(cl, c("ratioGCcontent"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
GCcontentByGene <- foreach(i=1:length(rownames(exprSQ)),.combine=rbind,.verbose=F)%dopar%ratioGCcontent(rownames(exprSQ)[i],GCcontentTab)
end <- Sys.time()
stopCluster(cl)
cat(paste("GC content ratios calculated in ",end-start,"\n"))
rm(cl,end,start,GCcontentTab,exonsdef)


GCcontent <- GCcontentByGene
rownames(GCcontent) <- GCcontent[,1]
GCcontent <- as.data.frame(GCcontent[,-1]) 
length <- length[as.character(rownames(GCcontent))]  
stopifnot(identical(names(length),rownames(GCcontent)))
GCcontent <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
colnames(GCcontent) <- c("GCcontent","length")
# we update gene expression with the filtered genes
rm(length,geneLength,genesList)
GCcontent$GCcontent <- as.numeric(GCcontent$GCcontent)

library(cqn)
library(scales)

## CQN
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),rownames(GCcontent)))              

cat("Conditional quantile normalisation \n")
my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent,sizeFactors=librarySize, verbose = TRUE)




GCcontent <- read.delim("data/general/GCcontRegionsExonic",row.names=1,sep=" ")
# remove the NaN in the GC content matrix
GCcontent[GCcontent=="NaN"]=""
colnames(GCcontent) <- "GCcontent"
GCcontent$GCcontent <- as.numeric(GCcontent$GCcontent)*100





