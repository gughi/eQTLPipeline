## main for genic quantification only exonic

rm(list=ls())
setwd("/home/guelfi/eQTLPipeline")
sink("logExons.log")
library(devtools)
load_all()

## Now we correct for PEER using simple quantification exons 
load("data/expr/rawCounts/genic/exons.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")

## dim(countsTable)
## [1] 639478    180

## convert the genes that have NAs
countsTable[is.na(countsTable)]=0
## remove genes that not expressed in any gene
countsTable <- countsTable[rowSums(countsTable>0)>0,]

## dim(countsTable)
## [1] 564270    180

cat("Processing PUTM \n")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]

# now we select the expression for the PUTM only samples+
expr <- countsTable[,as.character(PUTM$A.CEL_file)]
rm(countsTable)

## dim(expr)
## [1] 564270    105

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content genic and gene length

#geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
#geneLength <- geneLength[as.character(rownames(expr)),c(3,1:2)]

# load the definition of the exons
transcriptomeInfo <- read.csv("data/general/transcriptomeInfo.csv")
rownames(transcriptomeInfo) <- transcriptomeInfo$names
length <- transcriptomeInfo[as.character(rownames(expr)),3]
names(length) <-  as.character(rownames(expr))
## calculation of genes only exons length
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

save(RPKM.std,file="data/expr/rawCounts/genic/exons.RPKM.PUTM.rda")

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std)

## filtering
## write log
cat(paste("Number of Genes after filtering:",length(genesList),"\n"))
rm(RPKM.std)
expr <- expr[as.character(genesList),]

## obtain the chromosome number
genes <- unlist(lapply(strsplit(as.character(rownames(expr)),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

library("biomaRt")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


geneNames <- getBM(attributes=c("ensembl_gene_id","chromosome_name"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=genes, mart=ensembl)

rownames(geneNames) <- geneNames$ensembl_gene_id
geneNames$ensembl_gene_id <- NULL

exondef <- geneNames[as.character(genes),]
exondef <- cbind(row.names(expr),exondef)
exondef <- cbind(exondef,transcriptomeInfo[as.character(exondef[,1]),])
exondef <- exondef[,c(2:4,6)]
colnames(exondef) <- c("Chromosome.Name","Exon.Chr.Start..bp.","Exon.Chr.End..bp.","Ensembl.Gene.ID")



exons <- apply(exondef,1,function(x){paste(x[1],as.integer(x[2]),as.integer(x[3]),x[4],sep='\t')})


write.table(data.frame(exons), file = paste0("data/general/exons.BED"), row.names = F, 
            col.names = F, quote = F)

## we filter things that not match with the fasta file

##system("grep -v HG* data/general/exons.BED  | grep -v LRG* | grep -v HS* | cat > data/general/exons.BED ")

cmd <- paste0("/apps/BEDTools/2.24.0/bin/bedtools nuc -fi /home/ukbec/bowtie2Index/genome37.72.fa -bed data/general/exons.BED > data/general/GCexons")

## calculate GC content with bedtools
system(cmd)

GCcontentTab <- read.delim("data/general/GCexons")
cat(paste("GC content saved in data/general/GCexons","\n"))

rm(cmd,genes,ensembl,genesList)
## detectCores()
## [1] 24
notPresent <- rownames(expr)[!rownames(expr) %in% GCcontentTab$X4_usercol]
Sys.time()
GCcontentByExons <- sapply(notPresent, function(x)
    {
      if(x %in% GCcontentTab$X4_usercol)
      {  
        return(GCcontentTab[which(GCcontentTab$X4_usercol %in% x),"X6_pct_gc"])
      }else{
        tmp <- unlist(strsplit(x,":"))
        
        GCtmp <- GCcontentTab[grep(as.character(tmp[1]),GCcontentTab$X4_usercol,fixed = TRUE),]
        
        if (nrow(GCtmp)>0)
        {
          if(which.min(c(min(abs(GCtmp$X2_usercol-exondef[x,2])),min(abs(GCtmp$X3_usercol-exondef[x,2]))))==1)
          {
            return(GCtmp[which.min(abs(GCtmp$X2_usercol-exondef[x,2])),"X6_pct_gc"])
          
          }else{
            return(GCtmp[which.min(abs(GCtmp$X3_usercol-exondef[x,2])),"X6_pct_gc"])
        
          }
        }else{
          
          ## remove exons in case don't have any exon for that gene
          ## with -1 mean we can't calculate GC ocntent
          return(-1)
          
        }
      }
      
    }
    )

toRemove <- GCcontentByExons[which(GCcontentByExons < 0) ]
GCcontentByExons <- GCcontentByExons[-which(names(GCcontentByExons) %in% names(toRemove))]
## remove from the expression the exons that we eren't able to calculate the GC rate( There were 6 exons)
expr <- expr[-which(rownames(expr) %in% names(toRemove)),]

GCcontent <- as.data.frame(cbind(as.character(GCcontentTab$X4_usercol),GCcontentTab$X6_pct_gc))
rownames(GCcontent) <-GCcontent$V1
GCcontent <- GCcontent[as.character(rownames(expr)),]
GCcontent$V1 <- NULL
colnames(GCcontent)<- "GCcontent"
GCcontentByExons <- as.data.frame(GCcontentByExons)
colnames(GCcontentByExons) <- "GCcontent"
GCcontent <- rbind(GCcontent,GCcontentByExons)
GCcontent <- as.data.frame(GCcontent[as.character(rownames(expr)),])
rownames(GCcontent)<- rownames(expr)
rm(GCcontentByExons)


##rownames(GCcontent) <- GCcontent[,1]
##GCcontent <- as.data.frame(GCcontent[,-1]) 
length <- length[as.character(rownames(GCcontent))]  
stopifnot(identical(names(length),rownames(GCcontent)))

GCcontent <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
colnames(GCcontent) <- c("GCcontent","length")
# we update gene expression with the filtered genes
rm(length,geneLength,genesList,toRemove)
GCcontent$GCcontent <- as.numeric(GCcontent$GCcontent)

library(cqn)
library(scales)

## CQN
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),rownames(GCcontent)))              

cat("Conditional quantile normalisation \n")
my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent,sizeFactors=librarySize, verbose = TRUE)

png(paste0("plots/exons/CQNSNIG.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset


cat(paste("Number of Genes and samples",dim(RPKM.cqn),"\n"))
# 275703 exons in the analysis

# save results
cat("Saving the the RPKM CQN normalised in data/expr/normalisedCounts/genic/exons/RPKM.cqn.PUTM")
##save(RPKM.cqn,file="data/expr/normalisedCounts/exonsRPKM.cqn.rda",compress="bzip2")

PUTM$U.Region_simplified <- NULL
covs <- PUTM 
rownames(covs) <- covs$A.CEL_file
#convert the female and male info in numeric
covs[covs=="M"]=0
covs[covs=="F"]=1
covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
covs[,5] <- as.numeric(covs[,5])
covs[,6] <- as.numeric(covs[,6])
lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
rownames(lanes) <- gsub("CEL","",rownames(lanes))
covs <- cbind(covs,librarySize[as.character(rownames(covs))])
covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                    "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")

save(RPKM.cqn,PUTM,covs,file="data/expr/normalisedCounts/genic/exons/RPKM.cqn.PUTM")


### Residual correction ###

rm(list=ls())

load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.PUTM")

##doSwamp(RPKM.cqn,covs)

PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                               "data/expr/normalisedCounts/genic/exons/resids.PUTM.rda")

##doSwamp(resids,covs)

############
### SNIG ###    
############


rm(list=ls())
setwd("/home/guelfi/eQTLPipeline")
library(devtools)
load_all()

## Now we correct for PEER using simple quantification exons 
load("data/expr/rawCounts/genic/exons.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")

## dim(countsTable)
## [1] 639478    180

## convert the genes that have NAs
countsTable[is.na(countsTable)]=0
## remove genes that not expressed in any gene
countsTable <- countsTable[rowSums(countsTable>0)>0,]

## dim(countsTable)
## [1] 564270    180

cat("Processing SNIG \n")
SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),]

# now we select the expression for the SNIG only samples+
expr <- countsTable[,as.character(SNIG$A.CEL_file)]
rm(countsTable)

## dim(expr)
## [1] 564270    105

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(SNIG$A.CEL_file),]
names(librarySize) <- as.character(SNIG$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content genic and gene length

#geneLength <- read.delim("data/general/ensemblGenes.txt",row.names=1)
#geneLength <- geneLength[as.character(rownames(expr)),c(3,1:2)]

# load the definition of the exons
transcriptomeInfo <- read.csv("data/general/transcriptomeInfo.csv")
rownames(transcriptomeInfo) <- transcriptomeInfo$names
length <- transcriptomeInfo[as.character(rownames(expr)),3]
names(length) <-  as.character(rownames(expr))
## calculation of genes only exons length
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

save(RPKM.std,file="data/expr/rawCounts/genic/exons.RPKM.SNIG.rda")

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std)

## filtering
## write log
cat(paste("Number of Genes after filtering:",length(genesList),"\n"))
rm(RPKM.std)
expr <- expr[as.character(genesList),]

## obtain the chromosome number
genes <- unlist(lapply(strsplit(as.character(rownames(expr)),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

library("biomaRt")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


geneNames <- getBM(attributes=c("ensembl_gene_id","chromosome_name"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=genes, mart=ensembl)

rownames(geneNames) <- geneNames$ensembl_gene_id
geneNames$ensembl_gene_id <- NULL

exondef <- geneNames[as.character(genes),]
exondef <- cbind(row.names(expr),exondef)
exondef <- cbind(exondef,transcriptomeInfo[as.character(exondef[,1]),])
exondef <- exondef[,c(2:4,6)]
colnames(exondef) <- c("Chromosome.Name","Exon.Chr.Start..bp.","Exon.Chr.End..bp.","Ensembl.Gene.ID")



exons <- apply(exondef,1,function(x){paste(x[1],as.integer(x[2]),as.integer(x[3]),x[4],sep='\t')})


write.table(data.frame(exons), file = paste0("data/general/exons.BED"), row.names = F, 
            col.names = F, quote = F)

## we filter things that not match with the fasta file

##system("grep -v HG* data/general/exons.BED  | grep -v LRG* | grep -v HS* | cat > data/general/exons.BED ")

cmd <- paste0("/apps/BEDTools/2.24.0/bin/bedtools nuc -fi /home/ukbec/bowtie2Index/genome37.72.fa -bed data/general/exons.BED > data/general/GCexons")

## calculate GC content with bedtools
system(cmd)

GCcontentTab <- read.delim("data/general/GCexons")
cat(paste("GC content saved in data/general/GCexons","\n"))

rm(cmd,genes,ensembl,genesList)
## detectCores()
## [1] 24
notPresent <- rownames(expr)[!rownames(expr) %in% GCcontentTab$X4_usercol]
Sys.time()
GCcontentByExons <- sapply(notPresent, function(x)
{
  if(x %in% GCcontentTab$X4_usercol)
  {  
    return(GCcontentTab[which(GCcontentTab$X4_usercol %in% x),"X6_pct_gc"])
  }else{
    tmp <- unlist(strsplit(x,":"))
    
    GCtmp <- GCcontentTab[grep(as.character(tmp[1]),GCcontentTab$X4_usercol,fixed = TRUE),]
    
    if (nrow(GCtmp)>0)
    {
      if(which.min(c(min(abs(GCtmp$X2_usercol-exondef[x,2])),min(abs(GCtmp$X3_usercol-exondef[x,2]))))==1)
      {
        return(GCtmp[which.min(abs(GCtmp$X2_usercol-exondef[x,2])),"X6_pct_gc"])
        
      }else{
        return(GCtmp[which.min(abs(GCtmp$X3_usercol-exondef[x,2])),"X6_pct_gc"])
        
      }
    }else{
      
      ## remove exons in case don't have any exon for that gene
      ## with -1 mean we can't calculate GC ocntent
      return(-1)
      
    }
  }
  
}
)

toRemove <- GCcontentByExons[which(GCcontentByExons < 0) ]
GCcontentByExons <- GCcontentByExons[-which(names(GCcontentByExons) %in% names(toRemove))]
## remove from the expression the exons that we eren't able to calculate the GC rate( There were 6 exons)
expr <- expr[-which(rownames(expr) %in% names(toRemove)),]

GCcontent <- as.data.frame(cbind(as.character(GCcontentTab$X4_usercol),GCcontentTab$X6_pct_gc))
rownames(GCcontent) <-GCcontent$V1
GCcontent$V1 <- NULL
colnames(GCcontent)<- "GCcontent"
GCcontentByExons <- as.data.frame(GCcontentByExons)
colnames(GCcontentByExons) <- "GCcontent"
GCcontent <- rbind(GCcontent,GCcontentByExons)
GCcontent <- as.data.frame(GCcontent[as.character(rownames(expr)),])
rownames(GCcontent)<- rownames(expr)
rm(GCcontentByExons)


##rownames(GCcontent) <- GCcontent[,1]
##GCcontent <- as.data.frame(GCcontent[,-1]) 
length <- length[as.character(rownames(GCcontent))]  
stopifnot(identical(names(length),rownames(GCcontent)))

GCcontent <- cbind(GCcontent,length[as.character(rownames(GCcontent))])    
colnames(GCcontent) <- c("GCcontent","length")
# we update gene expression with the filtered genes
rm(length,geneLength,genesList,toRemove)
GCcontent$GCcontent <- as.numeric(GCcontent$GCcontent)

library(cqn)
library(scales)

## CQN
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),rownames(GCcontent)))              

cat("Conditional quantile normalisation \n")
my.cqn <- cqn(expr, lengths = GCcontent$length,x = GCcontent$GCcontent,sizeFactors=librarySize, verbose = TRUE)

png(paste0("plots/exons/CQNSNIG.jpeg"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset


cat(paste("Number of Genes and samples",dim(RPKM.cqn),"\n"))
# 275703 exons in the analysis

# save results
cat("Saving the the RPKM CQN normalised in data/expr/normalisedCounts/genic/exons/RPKM.cqn.SNIG")
##save(RPKM.cqn,file="data/expr/normalisedCounts/exonsRPKM.cqn.rda",compress="bzip2")

SNIG$U.Region_simplified <- NULL
covs <- SNIG 
rownames(covs) <- covs$A.CEL_file
#convert the female and male info in numeric
covs[covs=="M"]=0
covs[covs=="F"]=1
covs <- as.data.frame(apply(covs[,c(2:5,7:9)], 2, as.factor))
covs[,c(1:4,7)] <- as.data.frame(apply(covs[,c(1:4,7)], 2, as.numeric))
covs[,5] <- as.numeric(covs[,5])
covs[,6] <- as.numeric(covs[,6])
lanes <- read.csv("data/general/QCmetrics.csv",row.names=8)
rownames(lanes) <- gsub("CEL","",rownames(lanes))
covs <- cbind(covs,librarySize[as.character(rownames(covs))])
covs <- cbind(covs,lanes[as.character(rownames(covs)),c(9,19,20,25)])
colnames(covs) <- c("Age","PMI","RIN","Gender","CODE","OVation_Batch",
                    "TotReadsNoAdapt","LibrarySize","LanesBatch","uniqueMappedRead","FragLengthMean","ExonicRate")

save(RPKM.cqn,SNIG,covs,file="data/expr/normalisedCounts/genic/exons/RPKM.cqn.SNIG")



### PCA after normalisation ###

load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.SNIG")
RPKM.cqn.SNIG <- RPKM.cqn
load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.PUTM")
RPKM.cqn.PUTM <- RPKM.cqn

rm(RPKM.cqn)

comJunc <- intersect(rownames(RPKM.cqn.SNIG),rownames(RPKM.cqn.PUTM))
length(comJunc)
RPKM.cqn <- cbind(RPKM.cqn.SNIG[as.character(comJunc),],RPKM.cqn.PUTM[as.character(comJunc),])

PCAres<- prcomp(t(RPKM.cqn))
par(mfrow=c(1,1))

plot(PCAres, main="PCA axes exons(PUTM + SNIG)")
plot(PCAres$x[,1],PCAres$x[,2],main = "PC1 vs PC2 by tissue (exons)",xlab="PC1",ylab="PC2", ylim=c(-350,250),xlim=c(-500,850))

points(PCAres$x[PUTM$A.CEL_file,1],PCAres$x[PUTM$A.CEL_file,2],col="red")
points(PCAres$x[SNIG$A.CEL_file,1],PCAres$x[SNIG$A.CEL_file,2],col="blue")
legend("bottomright", c("PUTM", "SNIG"), pch = 1,col=c("red","blue"),title="tissue")    



### Residual correction ###

rm(list=ls())

load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.SNIG")

##doSwamp(RPKM.cqn,covs)

PEERRNDPEER18 <- read.csv("testPEER/RNDMPEER18",row.names=1)
PEERRNDPEER18 <- PEERRNDPEER18[colnames(RPKM.cqn),c(1:2,4:16)]

resids <- doResidualCorrection(t(RPKM.cqn),PEERRNDPEER18,
                               "data/expr/normalisedCounts/genic/exons/resids.SNIG.rda")


setwd("/home/guelfi/eQTLPipeline")

### Run the eQTL analysis
writeSH(nameSH="runCisEQTL.sh",logName="runCisEQTL",
        cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/parRunCiseQTLExons.R"),numThreads=8)


system("qsub runCisEQTL.sh")

 
### Run the Sentinalisation
writeSH(nameSH="LDsentinalisation.sh",logName="LDsentinalisation",
        cmdScript=paste0("/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=",getwd(),"/sentiExons.R"),numThreads=8)

system("qsub LDsentinalisation.sh")

# 
# ## collect all the eQTL in one file
# ## Collect for PUTM
# system("perl getAllEQTLsent.pl data/results/genic/geneExons/resMatrixEQTL/sentinalised/PUTM data/results/finaleQTLs/geneExonic.PUTM.txt")
# ## Collect for SNIG
# system("perl getAllEQTLsent.pl data/results/genic/geneExons/resMatrixEQTL/sentinalised/SNIG data/results/finaleQTLs/geneExonic.SNIG.txt")
# 
# # number of eTQL in PUTM
# system("wc -l data/results/finaleQTLs/geneExonic.PUTM.txt")
# # 1609
# 
# # number of eTQL in SNIG
# system("wc -l data/results/finaleQTLs/geneExonic.SNIG.txt | wc -l")
# # 951
# 
# 
# eQTLPUTM <- read.delim("data/results/finaleQTLs/geneExonic.PUTM.txt",sep=" ")
# table(eQTLPUTM$myFDR<0.05)
# #     FALSE  TRUE
# #     373  1235
# table(eQTLPUTM$myFDR<0.01)
# #     FALSE  TRUE
# #     841   767
# 
# 
# eQTLSNIG <- read.delim("data/results/finaleQTLs/geneExonic.SNIG.txt",sep=" ")
# table(eQTLSNIG$myFDR<0.05)
# #     FALSE  TRUE
# #     300   650
# table(eQTLSNIG$myFDR<0.01)
# #     FALSE  TRUE
# #     618   332
# 
# ###########################
# ####    Annotation ########
# ###########################
# 
# ### PUTM ###
# 
# 
# library("biomaRt")
# ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
#                    dataset="hsapiens_gene_ensembl")
# 
# 
# geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","gene_biotype","description"),
#                    verbose = T,
#                    filters="ensembl_gene_id",
#                    values=eQTLPUTM$gene, mart=ensembl)
# 
# 
# 
# rownames(geneNames)<- geneNames$ensembl_gene_id
# geneNames$ensembl_gene_id <- NULL
# 
# 
# eQTLPUTM <- cbind(eQTLPUTM,geneNames[as.character(eQTLPUTM$gene),])
# 
# save(eQTLPUTM,file="data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
# write.csv(eQTLPUTM,file="data/results/finaleQTLs/geneExonic.Ann.PUTM.csv",row.names=F)    
# 
# ### SNIG ###
# 
# rm(geneNames)
# 
# geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","start_position","end_position","gene_biotype","description"),
#                    verbose = T,
#                    filters="ensembl_gene_id",
#                    values=eQTLSNIG$gene, mart=ensembl)
# 
# 
# 
# rownames(geneNames)<- geneNames$ensembl_gene_id
# geneNames$ensembl_gene_id <- NULL
# 
# 
# eQTLSNIG <- cbind(eQTLSNIG,geneNames[as.character(eQTLSNIG$gene),])
# 
# save(eQTLSNIG,file="data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
# write.csv(eQTLSNIG,file="data/results/finaleQTLs/geneExonic.Ann.SNIG.csv",row.names=F)    
# 
# 
# 
# 




