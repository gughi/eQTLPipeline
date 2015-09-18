

load("data/expr/rawCounts/genic/exprSQ.rda")
exprSQ[is.na(exprSQ)]=0
## remove genes that not expressed in any gene
exprSQ <- exprSQ[rowSums(exprSQ>0)>0,]

load("data/general/sampleInfo.rda")



PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),]

load("data/general/genesWidthExonic.rda")

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[colnames(exprSQ),]
names(librarySize) <- as.character(colnames(exprSQ))

expr <- exprSQ
length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[as.character(rownames(expr))]
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

library(easyRNASeq)
RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]


ensemblRef <- read.delim(file="data/general/ensemblGenes.txt", as.is=T,header=T)
#order by position
ensemblRef <- ensemblRef[which(ensemblRef$Ensembl.Gene.ID %in% rownames(RPKM.std)),]
geneList <- ensemblRef[order(ensemblRef$Chromosome.Name,ensemblRef$Gene.Start..bp.),1]
RPKM.std <- RPKM.std[as.character(geneList),as.character(c(PUTM$A.CEL_file,SNIG$A.CEL_file))]

library(dplyr)
library(NMF)
library(RColorBrewer)

aheatmap(t(RPKM.std[1:100,]), color = "-RdBu:50", scale = "col", breaks = 0,
         annRow = c(PUTM$U.Region_simplified,SNIG$U.Region_simplified), annColors = "Set2", 
         treeheight=c(200, 50),distfun = "pearson",
         fontsize=13, cexCol=.7, 
         filename="heatmap.png", width=8, height=16)

heatmap(t(RPKM.std[1:100,]), color = "-RdBu:50", scale = "col", breaks = 0,
         annRow = c(PUTM$U.Region_simplified,SNIG$U.Region_simplified), annColors = "Set2", 
         treeheight=c(200, 50),
         fontsize=13, cexCol=.7, 
         filename="heatmap.png", width=8, height=16)

heatmap(t(RPKM.std[1:10,]))
install.packages("gplots")
library("gplots")

heatmap.2(t(RPKM.std[1:110,]),trace = "none")

head(RPKM.std)

iris2 = iris # prep iris data for plotting
rownames(iris2) = make.names(iris2$Species, unique = T)
iris2 = iris2 %>% select(-Species) %>% as.matrix()
heatmap(iris2, color = "-RdBu:50", scale = "col", breaks = 0,
         annRow = iris["Species"], annColors = "Set2", 
         distfun = "pearson", treeheight=c(200, 50), 
         fontsize=13, cexCol=.7, 
         filename="heatmap.png", width=8, height=16)


