## collect the exon reads in bioscience


library("DEXSeq")
tmpList <- read.csv("/home/skgtmgu/DEXSeqFiles.csv",header=F)
tmpList <-  gsub("/scratch2/","/cluster/scratch3/",tmpList[,1])   ## this is because they mounted the disk in the new cluster
flattened <- "/home/skgtmgu/transcriptome/Homo_sapiens.ExonLevel.GRCh37.72.gff"
sampleInfo <- read.csv("/home/skgtmgu/samplesInfo.csv",header=T,row.names=7)

sampleID <- strsplit(as.character(tmpList),"/",fixed=TRUE)
sampleID <- unlist(lapply(sampleID, function(x){x[6]}))
tmpList <- as.data.frame(tmpList)
tmpList$V2 <- sampleID

tmpList$V2 <- gsub("Sample_","",tmpList$V2)
tmpList$V2 <- gsub("CEL","",tmpList$V2)
rm(sampleID)
sampleTable <-data.frame(cbind(tmpList$V2,as.character(sampleInfo[tmpList$V2,6])))

rownames(sampleTable) <- sampleTable$X1
sampleTable$X1 <- NULL
colnames(sampleTable) <- "tissue"
dxd <- DEXSeqDataSetFromHTSeq(as.character(tmpList$tmpList),sampleData=sampleTable,design= ~ sample + exon + tissue:exon,flattenedfile=flattened)




## merge bioscience and apollo counts
countsTableApollo <- read.csv("/home/guelfi/expressionData/exonQuantification/apolloExpr.csv",row.names=1)
countsTableBio <- read.csv("/home/guelfi/expressionData/exonQuantification/bioExpr.csv",row.names=1)

countsTable <- cbind(countsTableApollo,countsTableBio)
rm(countsTableApollo,countsTableBio)


save("countsTable",file="data/expr/rawCounts/genic/exons.rda")







