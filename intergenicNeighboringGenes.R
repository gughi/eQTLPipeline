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
                                                                                                                 names = eQTLPUTM[i,"gene"])),GR)
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




## now we look correlation to exons

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

## we definde the transcriptome wiht exons

flattenedfile <- "data/general/Homo_sapiens.ExonLevel.GRCh37.72.gff"
aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, header = FALSE)
colnames(aggregates) <- c("chr", "source", "class", "start", 
                          "end", "ex", "strand", "ex2", "attr")
aggregates$strand <- gsub("\\.", "*", aggregates$strand)
aggregates <- aggregates[which(aggregates$class == "exonic_part"),] 
aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
transcripts <- strsplit(transcripts, "\\+")
exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1",aggregates$attr)
exoninfo <- GRanges(as.character(aggregates$chr), 
                    IRanges(start = aggregates$start, 
                            end = aggregates$end), strand = aggregates$strand)
names(exoninfo) <- paste(aggregates$gene_id, exonids,sep = ":E")
rm(flattenedfile,exonids)

distanceNearest <- function(GRQuery,GR){
  
  x <- countOverlaps(GRQuery,GR,ignore.strand=T)
  idx <- nearest(GRQuery,GR,ignore.strand=T)
  dis <- distanceToNearest(GRQuery,GR,ignore.strand=T)
  nearGene <- names(GR[idx,])
  stra <- as.character(strand(exoninfo[idx,]))
  return(c(distance=as.data.frame(dis)$distance,gene=nearGene,overlap=x,strand=stra))
}

library(doParallel)
library(foreach)

detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(4)
clusterExport(cl, c("overlapsAny","countOverlaps","distanceToNearest","GRanges","IRanges","Rle","nearest","strand"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
distance <- foreach(i=1:nrow(eQTLPUTM),.combine=rbind,.verbose=F)%dopar%distanceNearest(GRanges(seqnames = Rle(gsub("chr","",eQTLPUTM[i,"chr"])),
                                                                                                ranges = IRanges(start=eQTLPUTM[i,"start"],
                                                                                                                 end = eQTLPUTM[i,"end"],
                                                                                                                 names = eQTLPUTM[i,"gene"])),exoninfo)
end <- Sys.time()
end-start
stopCluster(cl)
table(distance[,1]>"0")
rm(cl,GR,end,start)

eQTLPUTM <- cbind(eQTLPUTM,distance)
colnames(eQTLPUTM)[2] <- "regionID"

load("data/expr/rawCounts/genic/exons.rda")
rm(RPKM.cqn.PUTM)
load("data/expr/rawCounts/intergenic/exprIntergenic.PUTM.rda")


load("data/general/sampleInfo.rda")
## we select only samples from PUTM
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]  

correla <- NULL
for (i in 1:nrow(eQTLPUTM))
{
  correla[i] <-  cor(as.numeric(exprIntergenic[as.character(eQTLPUTM$regionID[i]),as.character(PUTM$A.CEL_file)]),
                     as.numeric(countsTable[as.character(eQTLPUTM$gene[i]),as.character(PUTM$A.CEL_file)]))
  
}

eQTLPUTM <- cbind(eQTLPUTM,corre=correla)

eQTLPUTMTmp <- eQTLPUTM[order(eQTLPUTM$distance),]

## removing regions that have correlation equal to 0 to visualisation perposes 
eQTLPUTMTmp <- eQTLPUTMTmp[-which(is.na(eQTLPUTMTmp$corre)),]

plot(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre^2,las=2, main="correlation and distances (intergenic vs exons quantification)",
     xlab="distance",ylab="r^2")

hist(as.numeric(distance[which(as.numeric(distance[,1])>0),1]),breaks=40, main="histogram of distance btw intergenic and nearest exon", 
     xlab=" distance in bp")

hist(as.numeric(distance[which(as.numeric(distance[,1])>0 & as.numeric(distance[,1])<10000),1]),breaks=40, main="histogram of distance btw intergenic and nearest exon", 
     xlab=" distance in bp")


cor(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre)


plot(as.numeric(as.character(eQTLPUTMTmp$distance)),eQTLPUTMTmp$corre^2,las=2, main="correlation and distances (intergenic vs exons quantification) by starnd",
      xlab="distance",ylab="r^2")
points(as.character(eQTLPUTMTmp[which(as.character(eQTLPUTMTmp$strand) == "+"),"distance"]),
       (eQTLPUTMTmp[which(as.character(eQTLPUTMTmp$strand) =="+"),"corre"]^2),col="red")

points(as.character(eQTLPUTMTmp[which(as.character(eQTLPUTMTmp$strand) == "-"),"distance"]),
       (eQTLPUTMTmp[which(as.character(eQTLPUTMTmp$strand) =="-"),"corre"]^2),col="blue")

legend(x = "topright",legend =c("Forward","Reverse"), fill=c("red","blue"),)




head(eQTLPUTMTmp)

genes <- unlist(lapply(strsplit(as.character(eQTLPUTMTmp$gene),":"),function(x){x[1]}))
#genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

load("data/results/finaleQTLs/exons.PUTM.txt")
length(unique(eQTLPUTMTmp$gene))

head(eQTLPUTMExons)

eQTLPUTMExons <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")



library("gplots")
## look overlap from the exons
venn(list(PUTMInter=as.character(unique(eQTLPUTMTmp$gene)),
          PUTMEx=as.character(unique(eQTLPUTMExons$gene))))

genesFromExons <- unlist(lapply(strsplit(as.character(eQTLPUTMExons$gene),":"),function(x){x[1]}))

# look at the overlap from any exon in the gene
venn(list(PUTMInter=unique(genes),
          PUTMEx=unique(genesFromExons)))




genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x}))
genesFromExons <- unlist(sapply(strsplit(as.character(genesFromExons),"+",fixed=T),function(x){x}))

venn(list(PUTMInter=unique(genes),
          PUTMEx=unique(genesFromExons)))

load("data/results/finaleQTLs/exons.unsentinalised.rda")

## we look at overlap snp+gene and snp+intergenic region 
venn(list(PUTMInter=paste0(eQTLPUTMTmp$snps,eQTLPUTMTmp$gene),
          PUTMEx=paste0(my.eQTLs.PUTM[,1],my.eQTLs.PUTM[,2])))



## we now look at distance for all the intergenic regions

library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
rm(RPKM.cqn,covs,PUTM)

summary(starStopReg$width)

hist(starStopReg$width,main="length transcribed regions identified in PUTM")

library(GenomicRanges)

distanceNearest <- function(GRQuery,GR){
  
  x <- countOverlaps(GRQuery,GR,ignore.strand=T)
  idx <- nearest(GRQuery,GR,ignore.strand=T)
  dis <- distanceToNearest(GRQuery,GR,ignore.strand=T)
  nearGene <- names(GR[idx,])
  return(c(distance=as.data.frame(dis)$distance,gene=nearGene,overlap=x))
}


genesMap <- read.delim("data/general/ensemblGenes.txt")

head(genesMap)

## define the genomic regions
GR <- GRanges(seqnames = Rle(genesMap$Chromosome.Name),
              ranges = IRanges(start=genesMap$Gene.Start..bp.,
                               end = genesMap$Gene.End..bp.,
                               names = genesMap$Ensembl.Gene.ID))


library(doParallel)
library(foreach)


distanceNearest(GRanges(seqnames = Rle(gsub("chr","",as.character(starStopReg$chr[1]))),
                        ranges = IRanges(start=starStopReg$start[1],
                                         end =starStopReg$end[1]),
                        names=rownames(starStopReg)[1]),GR)


detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(6)
clusterExport(cl, c("overlapsAny","countOverlaps","distanceToNearest","GRanges","IRanges","Rle","nearest","strand"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
distance <- foreach(i=1:nrow(starStopReg),.combine=rbind,.verbose=F)%dopar%distanceNearest(GRanges(seqnames = Rle(gsub("chr","",as.character(starStopReg$chr[i]))),
                                                                                                ranges = IRanges(start=starStopReg$start[i],
                                                                                                                 end =starStopReg$end[i]),
                                                                                                names=rownames(starStopReg)[i]),GR)
end <- Sys.time()
end-start
stopCluster(cl)
table(distance[,1]>"0")
rm(cl,GR,end,start)


regInformation <- cbind(starStopReg,distance)

head(regInformation)

save(regInformation,file="data/general/defIntergenic.rda")


## this is to calculate the differences


#counts how many b are overlapping in two genomics regions
overlappingBP <- function(GR1,GR2)
{
  tmp <- intersect(GR1,GR2,ignore.strand=TRUE)
  if(length(tmp)>0){
    res <- abs((width(tmp) - width(GR1)))
    rm(tmp) 
  }else{
    
    res <- width(GR1)
    rm(tmp)
  }
  
  return(res)

}

# overlappingBP(GRanges(seqnames = Rle(gsub("chr","",regInformation$chr[i])),
#                       ranges = IRanges(start=regInformation$start[i],
#                                        end = regInformation$end[i])),
#               GR[as.character(regInformation$gene[i]),])



## we calculate the unique portion of our regions with the closest gene
cl <- makeCluster(6)
clusterExport(cl, c("intersect","width","GRanges","IRanges","Rle","length","abs"))
registerDoParallel(cl)
getDoParWorkers()
start <- Sys.time()
uniquePortion <- foreach(i=1:nrow(regInformation),.combine=c,.verbose=F)%dopar%overlappingBP(GRanges(seqnames = Rle(gsub("chr","",regInformation$chr[i])),
                                                                                                    ranges = IRanges(start=regInformation$start[i],
                                                                                                                     end = regInformation$end[i])),
                                                                                            GR[as.character(regInformation$gene[i]),])
end <- Sys.time()
end-start
stopCluster(cl)


head(regInformation)

regInformation <- cbind(regInformation,uniquePortion)


save(regInformation,file="data/general/defIntergenic.rda")



load(file="data/general/defIntergenic.rda")



