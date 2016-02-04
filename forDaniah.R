


library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)$pathname
setwd(dirname(path))
rm(path)

load("data/expr/rawCounts/genic/fullExExJun.rda")


LRRk2 <- "ENSG00000188906"

grep(LRRk2, as.character(exExJunAnn$geneID))

mapLRRK2 <- mapExon[grep(LRRk2, as.character(mapExon$geneID)),]


exonIDsTmp <- mapLRRK2$exonID


idxTmp <- unique(which(as.character(expr$Exon1ID) %in% as.character(exonIDsTmp)),
                 which(as.character(expr$Exon2ID) %in% as.character(exonIDsTmp)))


head(expr[,1:6])


exprForDaniah <- expr[idxTmp,]

which(rownames(exExJunAnn) %in% paste(exprForDaniah$Exon1ID,exprForDaniah$Exon2ID,sep = "_"))

library(devtools) 
load_all()

annForDaniah <- annExExJun(paste(exprForDaniah$Exon1ID,exprForDaniah$Exon2ID,sep = "_")[1],mapExon = mapExon)
for(i in 2:length(paste(exprForDaniah$Exon1ID,exprForDaniah$Exon2ID,sep = "_"))) 
{
  annForDaniah <- rbind(annForDaniah,annExExJun(paste(exprForDaniah$Exon1ID,exprForDaniah$Exon2ID,sep = "_")[i],mapExon = mapExon))
}

write.csv(annForDaniah,"tmp/exExForDaniahMap.csv")



