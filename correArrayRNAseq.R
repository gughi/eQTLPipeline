
load("data/expr/rawCounts/genic/exprSQ.rda")  ### load the expression
load("data/general/sampleInfo.rda")
mapGenes   <- read.csv(file="data/general/ensemblMapGenes.txt", header=T)


RNASeqExprPUTM <- exprSQ[,as.character(sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),"A.CEL_file"])]
RNASeqExprSNIG <- exprSQ[,as.character(sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),"A.CEL_file"])]

## get only the Genes that are in the MapGenes table some genesID weren't included in that table
RNASeqExprPUTM <- RNASeqExprPUTM[is.element(rownames(RNASeqExprPUTM),mapGenes$Ensembl.Gene.ID),]
RNASeqExprSNIG <- RNASeqExprSNIG[is.element(rownames(RNASeqExprSNIG),mapGenes$Ensembl.Gene.ID),]

GenesExpressedPUTM <-  mapGenes[is.element(mapGenes$Ensembl.Gene.ID,rownames(RNASeqExprPUTM)),2] 
GenesExpressedSNIG <-  mapGenes[is.element(mapGenes$Ensembl.Gene.ID,rownames(RNASeqExprSNIG)),2]

RNASeqExprPUTM <- RNASeqExprPUTM[as.character(mapGenes[match(unique(GenesExpressedPUTM), mapGenes[,2] ),1]),]
RNASeqExprSNIG <- RNASeqExprSNIG[as.character(mapGenes[match(unique(GenesExpressedSNIG),mapGenes[,2]),1]),]

colnames(RNASeqExprPUTM) <- gsub("/","_",sampleInfo[  match(colnames(RNASeqExprPUTM),sampleInfo[, 7]),1])
colnames(RNASeqExprSNIG) <- gsub("/","_",sampleInfo[  match(colnames(RNASeqExprSNIG),sampleInfo[, 7]),1])

rownames(RNASeqExprPUTM) <- unique(GenesExpressedPUTM)
rownames(RNASeqExprSNIG) <- unique(GenesExpressedSNIG)

## remove that genes that are not expressed in all the samples
RNASeqExprPUTM <- RNASeqExprPUTM[which(rowSums(RNASeqExprPUTM)>0),]
RNASeqExprSNIG <- RNASeqExprSNIG[which(rowSums(RNASeqExprSNIG)>0),]

cmd <- paste("cut -f1,6,7 data/general/t.map_Aug2012.txt")
mapGenestID <- read.table( pipe(cmd), header=T )
## Need to get the micro array gene's name in a matrix to then get only the first one.
microGenes <- t(sapply(strsplit(as.character(mapGenestID$Gene),","), '[', 1:max(sapply(strsplit(as.character(mapGenestID$Gene),","), length))))
length(which(unique(mapGenes$Associated.Gene.Name) %in% unique(microGenes[,1]))) ### Genes in common 11544
## get the genes that present in both databases
ensGenes <- mapGenes[which(unique(mapGenes$Associated.Gene.Name) %in% unique(microGenes[,1])),1] 
rm(cmd)

head(mapGenestID)

#################################
## get genes in common PUTM  ####
#################################


genesInCommon <-  microGenes[which(toupper(microGenes[,1]) %in% unique(GenesExpressedPUTM)),1]
## get only the unique genes
genesInCommon <- unique(genesInCommon)
tIDexpressed <-  mapGenestID[which(toupper(microGenes[,1]) %in% unique(GenesExpressedPUTM)),1]
### check how many genes we have in common between the two databases.
rm(GenesExpressedPUTM)
length(tIDexpressed)


tmp <- c()
for (i in 1:length(genesInCommon))
{
  ## get the transcript ID
  cmd <- paste("grep -iw  -e",genesInCommon[i], "data/general/t.map_Aug2012.txt","| cut -f1,7")
  df  <- read.table( pipe(cmd), header=F )
  rm(cmd)
  print(i)
  ## get the expression for tID getting the one that have maximum number of probes, at the end this is going to be the expression related with the gene
  tmp <- rbind(tmp, c(paste0("t",df[which.max(df[,2]),1]),genesInCommon[i]))
}


rm(i)
## create tmp file to use with fgrep
tmpf <- tempfile()
write.table(data.frame(tmp[,1]), file=tmpf, row.names=F, col.names=F, quote=F )
cmd <- paste("fgrep -e -iw -f", tmpf ,"/home/adai/DATA/expr/expr_PUTM.txt")
geneExprMicro <- read.table( pipe(cmd), header=F )  
rownames(geneExprMicro) <- tmp[match(geneExprMicro[,1],as.character(tmp[,1])),2]
##put the header
colnames(geneExprMicro) <- as.matrix(read.table( pipe("head -1 /home/adai/DATA/expr/expr_PUTM.txt"), header=F ))[1,]  
unlink(tmpf)
rm(cmd,tmp,tmpf,df)
## calculate the pearson's r

pearsonR <- c()

for (i in 1:length(colnames(RNASeqExprPUTM)))  
{
  test <- geneExprMicro[,colnames(RNASeqExprPUTM)[i]] 
  names(test) <- rownames(geneExprMicro)
  m <- as.matrix(cbind(test,RNASeqExprPUTM[names(test),colnames(RNASeqExprPUTM)[i]],deparse.level = 0))  
  tmpToremove <- apply(m,2, function(x) which(x == 0))
  ##remove the title
  ##plot(m[-tmpToremove[[2]],], main=paste(colnames(RNASeqExprPUTM)[i]), xlab="microArray",ylab="RNASeq",axes = FALSE)
  ##plot(m[-tmpToremove[[2]],], xlab="microArray",ylab="RNASeq",axes = FALSE)
  m[is.na(m)]=0
  pearsonR <- rbind(pearsonR,cor(as.vector(m[,1]),as.vector(m[,2])))
  print(paste(i,"of",length(colnames(RNASeqExprPUTM))))
  write.csv(pearsonR,file=paste0("data/general/pearson/",colnames(RNASeqExprPUTM)[i],"PUTM.csv"))
  rm(m)
}

hist(pearsonR,main="PUTM",breaks=50)
write.csv(pearsonR,file=paste0("data/general/pearson/pearsonRPUTM.csv"))






