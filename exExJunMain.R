
# Apollo
##setwd("/home/guelfi/eQTLPipeline")
# Caprica
setwd("/home/seb/projectsR/eQTLPipeline")
sink("logExonExonJunctions.log")
nCores <- 15
cat(paste("Number of cores",nCores,"\n"))
library(devtools)
load_all()

# Now we correct for PEER using simple quantification Exons+Introns



# expr <- read.delim("data/expr/rawCounts/genic/exonExonJunctions")
# ## first four columns have the exon1ID - exon2ID - chr and TSS
# 
# ## Loading information of exon -exon junctions
# # chr - start - end - ID
# map <- read.delim(pipe("grep GB data/general/exongrps.log | cut -f2-5 -d' '"),sep=" ",header=F)
# # Number of exons - length - max exon length If it's a complicated groupr and whether there are overlapping exons
# mapTmp <- read.delim(pipe("grep GI data/general/exongrps.log | cut -f2-6 -d' '"),sep=" ",header=F)
# map <- cbind(map,mapTmp)
# 
# ## load the map of the ID with the exon IDs
# mapTmp <- read.delim(pipe("grep GE data/general/exongrps.log"),sep=" ",header=F,fill=T , col.names=paste("col", 1:64, sep="_"))
# ##mapTmp[1,!is.na(mapTmp[1,])]
# 
# y <-1:nrow(mapTmp)
# IDs <- unlist(sapply(y, function(x){ rep(x,length(mapTmp[x,2:64][!is.na(mapTmp[x,2:64])]))}))
# 
# mapExon <- read.delim(pipe("grep -w E data/general/exongrps.log"),sep=" ",header=F)
# mapExon <- cbind(mapExon,IDs)
# ## columns check the Altrans manual if something unclear
# mapExon$V1 <- NULL
# 
# 
# ## clean the data
# 
# colnames(mapExon) <- c("chr","start", "end", "exonID","geneID", "transID", "unifiedExon",
#                      "differentStrand","strand","length","UnRegStartEnd",
#                      "groupID")
# 
# colnames(map) <- c("chr","start", "end", "groupID","noExons", "length",
#                    "maxExonLength","probleGroup","nonOverlaExon")
# 
# map$noExons <- gsub("NoExons:","",map$noExons)
# map$length <- gsub("Length:","",map$length)
# map$maxExonLength <- gsub("MaxExonLength:","",map$maxExonLength)
# save(map,mapExon,expr,file="data/expr/rawCounts/fullExExJun.rda")

load("data/expr/rawCounts/genic/fullExExJun.rda")

# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")
## reformat the sample names to match exon-exon Junction expression names

# investigate the PC1
# PCA PC1 vs PC2 
PCAres<- prcomp(t(expr[,5:ncol(expr)]))
par(mfrow=c(1,1))
PUTM <- sampleInfo[sampleInfo[, 6] == "PUTM",]
SNIG <- sampleInfo[sampleInfo[, 6] == "SNIG",]

plot(PCAres, main="PCA axis exon-exon juctions (PUTM + SNIG)")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue",xlab="PC1",ylab="PC2" )

exprSamNam <- sapply(PUTM$A.CEL_file,function(x){
  paste(unlist(strsplit(x,"_")),collapse="_")})

points(PCAres$x[paste0("Sample_",exprSamNam),1],
       PCAres$x[paste0("Sample_",exprSamNam),2],col="red")

exprSamNam <- sapply(SNIG$A.CEL_file,function(x){
  paste(unlist(strsplit(x,"_")),collapse="_")})

points(PCAres$x[paste0("Sample_",exprSamNam),1],
       PCAres$x[paste0("Sample_",exprSamNam),2],col="blue")

legend("bottomleft", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    

rm(PUTM,SNIG,PCAres)

doSwamp(t(expr[,5:ncol(expr)]),covs=sampleInfo)
## we calculate the GC content for all the junctions

write.table(data.frame(mapExon[,1:4]), file = paste0("data/general/regionJunction.BED"), row.names = F, 
            col.names = F, quote = F,sep="\t")

## we filter things that not match with the fasta file

system("grep -v HG* data/general/regionJunction.BED  | grep -v LRG* | grep -v HS* | cat > data/general/regionJunctionFiltered.BED")

cmd <- paste0("/apps/BEDTools/2.24.0/bin/bedtools nuc -fi /home/ukbec/bowtie2Index/genome37.72.fa -bed data/general/regionJunctionFiltered.BED > data/general/GCcontRegionJunctions; rm data/general/regionJunctionFiltered.BED ;rm data/general/regionJunction.BED")

## calculate GC content with bedtools
system(cmd)
rm(cmd,exprSamNam)

cat(paste("GC content saved in data/general/GCcontRegionJunctions","\n"))
GCcontentTab <- read.delim("data/general/GCcontRegionJunctions")



# now we select the expression for the PUTM only samples
exprAll <- expr
colnames(expr) <- paste0(gsub("Sample_","",colnames(expr)),"_")

# now we select the expression for the PUTM only samples
expr <- expr[,as.character(PUTM$A.CEL_file)]



## detectCores()
## [1] 24


librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content genic and gene length

juncdef <- exprAll[,1:4]

## calculation of genes only exons length

library(doParallel)
library(foreach)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(20)
clusterExport(cl, c("lengthJunction"))

registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
system.time(length <- foreach(i=1:nrow(juncdef[,1:2]),.combine=c,.verbose=F)%dopar%lengthJunction(juncdef[i,1:2],mapExon))
#geneswidth <- foreach(i=1:10,.combine=rbind,.verbose=F)%dopar%getRegionsWidth(rownames(expr)[i],exonsdef)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)


length <- apply(juncdef[,1:2],1,function(x) lengthJunction(x,mapExon))








 
  
  
  
  