# Differentially expressed analysis

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade") 
# biocLite("DESeq2")
# biocLite("easyRNASeq")


library(DESeq2)

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




GCcontent <- read.delim("ensemblRef/GCcontentExonic",row.names=1,sep=" ")
# remove the NaN in the GC content matrix
GCcontent[GCcontent=="NaN"]=""
colnames(GCcontent) <- "GCcontent"
GCcontent$GCcontent <- as.numeric(GCcontent$GCcontent)*100





