
library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

## load teh eQTL table
load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR <=0.05),1:7]

## exonJunc tells in how manysamples this regions has been seeing as an exon-junction
eQTLPUTM$exonJunc <- 0
## exonJuncCounts tells how many exon junctions has been seeying in total in all the samples
eQTLPUTM$exonJuncCounts <- 0

samples <- list.files("data/junctions/PUTM/")
samples <- unlist(lapply(strsplit(as.character(samples),"CEL_"),function(x){x[1]}))
samples <- paste0(unique(samples),"CEL")


## with this we get a summary for the novel exon in our eQTLs
getSummaryNovExon <- function(eQTL,intNovelExonPath)
{
  
  load(intNovelExonPath)
  
  intNovelExon$distance <- as.numeric(as.character(intNovelExon$distance))
  intNovelExon <- intNovelExon[which(intNovelExon$distance >0),]
  intNovelExon <- intNovelExon[which(as.character(intNovelExon$novelExon) == TRUE),]
  
  for(i in 1:nrow(eQTL))
  {
    if(is.element(as.character(eQTL$gene[i]),as.character(rownames(intNovelExon))))
    {
      eQTL$exonJunc[i] <- eQTL$exonJunc[i] +1
      eQTL$exonJuncCounts[i] <- eQTL$exonJuncCounts[i] + as.numeric(as.character(intNovelExon[as.character(eQTL$gene[i]),"counts"]))
    }
  }
  rm(intNovelExon)
  return(eQTL)
  
}


for(i in 1:length(samples))
{
  eQTLPUTM <- getSummaryNovExon(eQTLPUTM,intNovelExonPath=paste0("data/junctions/PUTM/",samples[i],"_intergenicNovelExon.rda"))
}

table(eQTLPUTM$exonJunc>0)
nrow(eQTLPUTM[which(eQTLPUTM$exonJunc>0),])

hist(eQTLPUTM[which(eQTLPUTM$exonJunc>0),"exonJunc"],breaks=40,main="Novel exon identify by sample",xlab="samples")

## example
regInformation["DER591",]

save(eQTLPUTM, file="data/results/finaleQTLs/intergenic.annNovelExon.PUTM.rda")

###########################
#### Beta interaction #####
###########################

load("data/general/defIntergenic.rda")

load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

eQTLs <- eQTLPUTM[which(eQTLPUTM$myFDR <=0.05),1:7]

regInformation <- regInformation[as.character(eQTLPUTM$gene),]

load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")
exprExonic <- resids
load("data/expr/normalisedCounts/intergenic/resids.PUTM.rda")
exprIntergenic <- resids
rm(resids,eQTLPUTM)

eQTLs <-  cbind(eQTLs,regInformation[as.character(eQTLs$gene),])
colnames(eQTLs)[2] <- "region"


## we now collect the information for each gene
eQTLExonic <- NULL
for(i in 1:nrow(eQTLs))
{
  if(is.element(as.character(eQTLs$gene[i]),colnames(exprExonic)))
  {
    tmp <- read.delim(paste0("data/results/genic/geneExons/fullResults/PUTM/",as.character(eQTLs$gene[i])),row.names = 1)
    eQTLExonic <- rbind(eQTLExonic,tmp[as.character(eQTLs$snps[i]),c("beta","t.stat","p.value","FDR")])
    rm(tmp)
    
  }else{
    eQTLExonic <- rbind(eQTLExonic,as.data.frame(t(as.matrix(c(beta=NA,t.stat=NA,p.value=NA,FDR=NA))),row.names=as.character(eQTLs$snps[i])))
  }
}
rm(i)


colnames(eQTLExonic) <- c("ge.beta","ge.t.stat","ge.p.value","ge.FDR")

eQTLs <- cbind(eQTLs,eQTLExonic)
rm(eQTLExonic)

library(devtools)
load_all()

f<- NULL
for (j in 1:nrow(eQTLs))
{
  ## load the expression for each individual gene for exonic and intronic
  if(is.element(as.character(eQTLs$gene[j]),colnames(exprExonic)))
  {
    exprE <- exprExonic[,as.character(eQTLs$gene[j])]
    exprInt <- exprIntergenic[,as.character(eQTLs$region[j])]
  
    ## add to the expression the tissue information or the exonic/intronic information
    expr <- cbind(exprE,"Exonic")
    exprtmp <- cbind(exprInt,"Intergenic")
    expr <- rbind(expr,exprtmp)
    rm(exprtmp,exprE,exprInt)
  
    ## load the snps
    load(paste0("/home/seb/eQTL/snps/byGene/",as.character(eQTLs$gene[j]),".rda"))
    
    ## load the samples information
    load("data/general/sampleInfo.rda")
    
    ## change the samples names and remove outlier
    sampleInfo$U.SD_No <- gsub("/","_",sampleInfo$U.SD_No)
    #   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "95_40"),]
    #   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "09_50"),]
    #   sampleInfo <- sampleInfo[-which(sampleInfo$U.SD_No %in% "034_07"),]
    #   
    rownames(sampleInfo) <- sampleInfo$A.CEL_file
    selSamples <- sampleInfo[as.character(rownames(expr)),1]
    
    ##snp0 <- markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
    
    snp0 <- markers[as.character(eQTLs$snps[j]),selSamples]
    ##snp0 <- round(snp0)
    
    ## markers[as.character(eQTLs$snps[j]),unique(as.character(gsub("/","_",sampleInfo$U.SD_No)))]
    
    ## rownames(sampleInfo) <- sampleInfo$A.CEL_file
    
    colnames(snp0) <- rownames(expr)
    f[j] <- betaInteraction(expr,snp0)
  }else
  {
    f[j] <- NA
  }
  
  rm(expr,markers,markers.info,snp0,selSamples,sampleInfo)
}

rm(j)
  
eQTLs <- cbind(eQTLs,p.val.bet.int=f)
rm(f)

## save(eQTLs,file="data/results/finaleQTLs/betaInterIntergenicNearestGene.rda")



load("data/results/finaleQTLs/betaInterIntergenicNearestGene.rda")
load(file="data/results/finaleQTLs/intergenic.annNovelExon.PUTM.rda")

eQTLs <- cbind(eQTLs,exonJunc=eQTLPUTM$exonJunc,exonJuncCounts=eQTLPUTM$exonJuncCounts)

head(eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),])


qq.plot(eQTLs[which(as.numeric(as.character(eQTLs$distance))<1),"p.val.bet.int"])
head(eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0),])


eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))==0 & as.numeric(as.character(eQTLs$distance)) >0),c("pvalue","ge.p.value","p.val.bet.int")]

eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0),"p.val.bet.int"]

table(is.na(eQTLs$p.val.bet.int))


qq.plot(eQTLs[which(as.numeric(as.character(eQTLs$ge.p.value))<1e-5),"p.val.bet.int"])


eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>50),]


## we calculate the fold change


load("data/expr/normalisedCounts/genic/geneExons/RPKM.cqn.PUTM")
exprExonic <- RPKM.cqn
load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
exprIntergenic <- RPKM.cqn
rm(RPKM.cqn,starStopReg,PUTM,covs)


for (j in 1:nrow(eQTLs))
{
  ## load the expression for each individual gene for exonic and intronic
  if(is.element(as.character(eQTLs$gene[j]),rownames(exprExonic)))
  {
    exprE <- exprExonic[as.character(eQTLs$gene[j]),]
    exprInt <- exprIntergenic[as.character(eQTLs$region[j]),]

    stopifnot(identical(names(exprInt),names(exprE)))
    
    
    RPKM.e <- mean(exprE)
    RPKM.i <- mean(exprInt)
    eQTLs[j,"fold.change"] <- RPKM.i/RPKM.e
    if(eQTLs[j,"fold.change"] < 1){
      eQTLs[j,"fold.change"] <- -1/eQTLs[j,"fold.change"]
    }
#     if(abs(eQTLs[j,"fold.change"]) > fold.theQTLshold & eQTLs[j,"FDRInter"] < p.val.theQTLshold){
#       eQTLs[j,"colors"] <- "red"  
#     }else
#       eQTLs[j,"colors"] <- "black"	
#       eQTLs$p.interaction
#     }
    
    eQTLs[j,"correlation"] <- cor(exprE,exprInt)
    rm(exprInt,exprE,RPKM.i,RPKM.e)
    }else
  {
    eQTLs[j,"fold.change"] <- NA
    eQTLs[j,"correlation"] <- NA
  }
  
 
}

head(eQTLs)

tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]

qq.plot(tmp$p.val.bet.int,main="Beta interaction test of potential novel exons")

head(tmp)

plot(log2(tmp$fold.change),log10(tmp$p.val.bet.int))
plot(tmp$correlation^2,-log10(tmp$p.val.bet.int))


plot(eQTLs$correlation^2,-log10(eQTLs$p.val.bet.int))

plot(log2(eQTLs$fold.change),-log10(eQTLs$p.val.bet.int))


tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>10),]


plot(log2(tmp$fold.change),-log10(tmp$p.val.bet.int))
plot(tmp$correlation^2,-log10(tmp$p.val.bet.int))

eQTLs[which(eQTLs$region == "DER591"),]



#### correlation and fold change with residuals
# 
# load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")
# exprExonic <- resids
# load("data/expr/normalisedCounts/intergenic/resids.PUTM.rda")
# exprIntergenic <- resids
# rm(resids)
# 
# 
# for (j in 1:nrow(eQTLs))
# {
#   ## load the expression for each individual gene for exonic and intronic
#   if(is.element(as.character(eQTLs$gene[j]),colnames(exprExonic)))
#   {
#     exprE <- exprExonic[,as.character(eQTLs$gene[j])]
#     exprInt <- exprIntergenic[,as.character(eQTLs$region[j])]
#     
#     stopifnot(identical(names(exprInt),names(exprE)))
#     
#     
#     RPKM.e <- mean(exprE)
#     RPKM.i <- mean(exprInt)
#     eQTLs[j,"fold.change"] <- RPKM.i/RPKM.e
#     if(eQTLs[j,"fold.change"] < 1){
#       eQTLs[j,"fold.change"] <- -1/eQTLs[j,"fold.change"]
#     }
#     #     if(abs(eQTLs[j,"fold.change"]) > fold.theQTLshold & eQTLs[j,"FDRInter"] < p.val.theQTLshold){
#     #       eQTLs[j,"colors"] <- "red"  
#     #     }else
#     #       eQTLs[j,"colors"] <- "black"	
#     #       eQTLs$p.interaction
#     #     }
#     
#     eQTLs[j,"correlation"] <- cor(exprE,exprInt)
#     rm(exprInt,exprE,RPKM.i,RPKM.e)
#   }else
#   {
#     eQTLs[j,"fold.change"] <- NA
#     eQTLs[j,"correlation"] <- NA
#   }
#   
#   
# }
# 
# 
# 
# tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
# plot(log2(tmp$fold.change),-log10(tmp$p.val.bet.int))
# plot(tmp$correlation^2,-log10(tmp$p.val.bet.int))
# 
# 
# plot(eQTLs$correlation^2,-log10(eQTLs$p.val.bet.int))
# 
# plot(log2(eQTLs$fold.change),-log10(eQTLs$p.val.bet.int))
# 


tmp <- eQTLs[which(as.numeric(as.character(eQTLs$uniquePortion))>0),]

tmp <- tmp[which(as.numeric(as.character(tmp$distance))==0),]

tmp <- tmp[which(as.numeric(as.character(tmp$overlap))<2),]


library(devtools)
setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()
qq.plot(tmp$p.val.bet.int,main="pval beta interaction intergenic regions attached to nearest gene")


plot(log2(tmp$fold.change),-log10(tmp$p.val.bet.int))
plot(tmp$correlation^2,-log10(tmp$p.val.bet.int))

regInformation["DER29903",]

tmp[which(tmp$p.val.bet.int == tmp$p.val.bet.int,na.rm = T),]

TMP2 <- tmp[which(tmp$p.val.bet.int > 0.05),]

mean(TMP2$uniquePortion/TMP2$width)

tmp3 <- tmp[which(tmp$p.val.bet.int < 0.05),]
t.test((tmp3$uniquePortion/tmp3$width),(TMP2$uniquePortion/TMP2$width))




#### After discussion we will separate the intergenic based on correlations and junctions.

load("data/results/finaleQTLs/betaInterIntergenicNearestGene.rda")
load(file="data/results/finaleQTLs/intergenic.annNovelExon.PUTM.rda")
eQTLs <- cbind(eQTLs,exonJunc=eQTLPUTM$exonJunc,exonJuncCounts=eQTLPUTM$exonJuncCounts)
rm(eQTLPUTM)

load("data/expr/rawCounts/intergenic/exprIntergenic.PUTM.rda")
exprIntergenic <- exprIntergenic[,5:ncol(exprIntergenic)]
load("data/expr/rawCounts/genic/exons.rda")
exprExons <- countsTable[,as.character(colnames(exprIntergenic))] 
rm(countsTable)


for (j in 1:nrow(eQTLs))
{
  ## load the expression for each individual gene for exonic and intronic
  exprE <- exprExons[grep(as.character(eQTLs$gene[j]),rownames(exprExons)),]
  exprInt <- exprIntergenic[as.character(eQTLs$region[j]),]
  stopifnot(identical(names(exprInt),names(exprE)))
    
  eQTLs[j,"correlation"] <- max(apply(exprE,1,function(x){cor(as.numeric(x),as.numeric(exprInt[1,]))^2} ))
  rm(exprInt,exprE,RPKM.i,RPKM.e)
  print(paste((j/nrow(eQTLs))*100 ,"% completed"))
  
}

## we filtered more the samples although this is not saved
eQTLs <- eQTLs[which(as.numeric(as.character(eQTLs$uniquePortion))>0),]


## this corresponds to scenario 1 intergenic linked to the neares gene with a junction
plot(as.numeric(as.character(eQTLs$distance)),as.numeric(as.character(eQTLs$correlation)), ylab = "Correlation",xlab="distance",main="Intergenic eQTLs - Distance Vs Correlation with nearest gene")
abline(h = 0.2,col="red")
abline(v=5000,col="red")

## scenario 2 misannotation
tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[!is.na(tmp$correlation),]
tmp <- tmp[which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]


points(as.numeric(as.character(tmp$distance)),as.numeric(as.character(tmp$correlation)),col="red")
rm(tmp)

## sceanrio3 indipendent intergenic
tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance))>0),]
tmp <- tmp[!is.na(tmp$correlation),]
##tmp <- tmp[which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- tmp[which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
points(as.numeric(as.character(tmp$distance)),as.numeric(as.character(tmp$correlation)),col="green")
rm(tmp)
## scenario 1
tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[!is.na(tmp$correlation),]
points(as.numeric(as.character(tmp$distance)),as.numeric(as.character(tmp$correlation)),col="blue")
rm(tmp)

legend("topright",legend = c("Novel exons","Misannotation","Independent intergenic","unknown"),col=c("blue","red","green","black"),pch = 1)
    
    

for (j in 1:nrow(eQTLs))
{
  ## load the expression for each individual gene for exonic and intronic
  exprE <- exprExons[grep(as.character(eQTLs$gene[j]),rownames(exprExons)),]
  exprInt <- exprIntergenic[as.character(eQTLs$region[j]),]
  stopifnot(identical(names(exprInt),names(exprE)))
  
  arrayTmp <- apply(exprE,1,function(x){cor(as.numeric(x),as.numeric(exprInt[1,]))^2})
  if(is.na(max(arrayTmp))){
    eQTLs[j,"exon"] <- NA
  }else{
    if(length(which(max(arrayTmp,na.rm = T)== arrayTmp))>1){
      eQTLs[j,"exon"] <- names(arrayTmp)[1]
    }else{
      eQTLs[j,"exon"] <- names(arrayTmp)[which(max(arrayTmp,na.rm = T) == arrayTmp)] 
    }
  }
  rm(exprInt,exprE,arrayTmp)
  print(paste((j/nrow(eQTLs))*100,"% completed"))
}
rm(j)

exprExonTmp <- exprExons[as.character(eQTLs$exon),]
exprExonTmp <- apply(exprExonTmp,1,mean)
exprExonTmp <- exprExons[!is.na(exprExonTmp[,1]),]
head(exprExonTmp)

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[!is.na(tmp$correlation),]
misann <- tmp[which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]
inIR <- tmp[which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
novelExons <- tmp[!is.na(tmp$correlation),]

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- tmp[!is.na(tmp$correlation),]


length(unique(inIR$exon))


toPlot <- rbind(exprExons[unique(as.character(misann$exon)),],exprExons[unique(as.character(inIR$exon)),]
                ,exprExons[unique(as.character(novelExons$exon)),],exprExons[unique(as.character(tmp$exon)),])


head(exprExons)
heatmap(as.matrix(toPlot))
library(gplots)
heatmap(as.matrix(toPlot),key=T,
          ##density.info="none", 
          ##trace="none", 
          RowSideColors = c(    # grouping row-variables into different
  rep("red", nrow(exprExons[unique(as.character(misann$exon)),])),   # categories, Measurement 1-3: green
  rep("green", nrow(exprExons[unique(as.character(inIR$exon)),])),    # Measurement 4-6: blue
  rep("blue", nrow(exprExons[unique(as.character(novelExons$exon)),])),
  rep("black", nrow(exprExons[unique(as.character(tmp$exon)),]))))


heatmap
 

heatmap.2(as.matrix(toPlot),
          col=redgreen(75), key=T, keysize=1.5,scale = "row",
          labRow=NA,ylab = "exons",symkey=F,trace="none",rowSideColors = c("red", "green","blue", "black")
          RowSideColors = c(    # grouping row-variables into different
            rep("red", nrow(exprExons[unique(as.character(misann$exon)),])),   # categories, Measurement 1-3: green
            rep("green", nrow(exprExons[unique(as.character(inIR$exon)),])),    # Measurement 4-6: blue
            rep("blue", nrow(exprExons[unique(as.character(novelExons$exon)),])),
            rep("black", nrow(exprExons[unique(as.character(tmp$exon)),]))))


par(lend = 1)

#Load necessary packages
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function

rowAnno = as.matrix(t(c(    # grouping row-variables into different
  rep("red", nrow(exprExons[unique(as.character(misann$exon)),])),   # categories, Measurement 1-3: green
  rep("green", nrow(exprExons[unique(as.character(inIR$exon)),])),    # Measurement 4-6: blue
  rep("blue", nrow(exprExons[unique(as.character(novelExons$exon)),])),
  rep("black", nrow(exprExons[unique(as.character(tmp$exon)),])))))

breaks=seq(-7, 7, by=0.1) #41 values
#now add outliers
breaks=append(breaks, 10)
breaks=append(breaks, -10, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="green",mid="black",high="red")

heatmap.2(as.matrix(toPlot),col=greenred(75), key=T, keysize=1.5,scale = "row",
          labRow=NA,ylab = "exons",trace="none",
            RowSideColors = rowAnno)

par(mar=c(0, 0, 0, 0))
heatmap.2(as.matrix(toPlot),col=mycol, key=T, keysize=1.5,scale = "",
            labRow=NA,ylab = "exons",symkey=F,trace="none",
            RowSideColors = rowAnno,breaks = breaks)


legend("topright",legend=c("misannotation", "ind. int. reg.", "novel Exon", "not defined"),
       fill=c("red", "green","blue", "black"), border=FALSE, bty="n", y.intersp = .7, cex=0.7)


breaks=seq(0, 10, by=1) #41 values
breaks=append(breaks, 10000)
breaks=append(breaks, -10, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="green",mid="black",high="red")


heatmap.2(as.matrix(toPlot),col=mycol, key=T, keysize=1.5,
          labRow=NA,ylab = "exons",symkey=F,trace="none", breaks=breaks,
          RowSideColors = rowAnno)

legend("topright",legend=c("misannotation", "ind. int. reg.", "novel Exon", "not defined"),
       fill=c("red", "green","blue", "black"), border=FALSE, bty="n", y.intersp = .7, cex=0.7)

## USING THE RPKM
load("data/expr/rawCounts/genic/exons.RPKM.PUTM.rda")

toPlot <- rbind(RPKM.std[unique(as.character(misann$exon)),],RPKM.std[unique(as.character(inIR$exon)),]
                ,RPKM.std[unique(as.character(novelExons$exon)),],RPKM.std[unique(as.character(tmp$exon)),])


heatmap.2(as.matrix(toPlot),
        col=redgreen(75), key=T, keysize=1.5,scale = "row",
        density.info="none", trace ="none",cexCol=0.9, labRow=NA,
        RowSideColors = c(    # grouping row-variables into different
          rep("red", nrow(RPKM.std[unique(as.character(misann$exon)),])),   # categories, Measurement 1-3: green
          rep("green", nrow(RPKM.std[unique(as.character(inIR$exon)),])),    # Measurement 4-6: blue
          rep("blue", nrow(RPKM.std[unique(as.character(novelExons$exon)),])),
          rep("black", nrow(RPKM.std[unique(as.character(tmp$exon)),]))))



## show the different distributions

head(exprExons)

hist(exprExons[unique(as.character(misann$exon)),])

exprExons[unique(as.character(inIR$exon)),]
exprExons[unique(as.character(novelExons$exon)),]
exprExons[unique(as.character(tmp$exon)),]



boxcox(c(0,0))


par(mar=c(3, 5, 5, 3))
boxplot(cbind(misAnnotation=log1p(as.numeric(unlist(exprExons[unique(as.character(misann$exon)),]))),
              indIntergenic=log1p(as.numeric(unlist(exprExons[unique(as.character(inIR$exon)),]))),
              novelExon=log1p(as.numeric(unlist(exprExons[unique(as.character(novelExons$exon)),]))),
              unknown=log1p(as.numeric(unlist(exprExons[unique(as.character(tmp$exon)),])))),col=c("red","green","blue","grey"),main="Raw counts",ylab="log(counts+1)")




## below we get the position of the region in respect of the nearest exon, in other words is it upstream or downstream



{

      
      library(biomaRt)
      ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                         dataset="hsapiens_gene_ensembl")
      
      geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","chromosome_name","start_position","end_position","strand"),
                         verbose = T,
                         filters="ensembl_gene_id",
                         values=eQTLs$gene, mart=ensembl)

      
      rownames(geneNames) <- geneNames$ensembl_gene_id  
      
      
      locReg <- NULL
      ## locReg contains the location of the intergenic in respect to the nearest gene
      
      for(i in 1:nrow(eQTLs))
      {
        gene <-geneNames[as.character(eQTLs$gene[i]),]
        if(gene$strand==1)
        {
          if(gene$end_position < eQTLs$end[i])
          {
            locReg[i] <- "down_stream"
          }else {
            locReg[i] <- "up_stream"
          }  
        }else{
          if(gene$start_position > eQTLs$start[i])
          {
            locReg[i] <- "down_stream"
          }else {
            locReg[i] <- "up_stream"
          }
        }
        
        rm(gene)
      }
  
}

eQTLs <- cbind(eQTLs,locReg)


save(eQTLs,file="data/results/novelIntergenicRegionsFinal/finalNovelIntergeniceQTL.PUTM.rda")



## we now plot upstream and down stream genes

plot(as.numeric(as.character(eQTLs$distance)),as.numeric(as.character(eQTLs$correlation)), ylab = "Correlation",xlab="distance",main="Distance vs Correlation by location from the nearest gene ")
abline(h = 0.2,col="red")
abline(v=5000,col="red")

## scenario 2 misannotation
points(as.numeric(as.character(eQTLs[which(as.character(eQTLs$locReg)=="down_stream"),"distance"])),as.numeric(as.character(eQTLs[which(as.character(eQTLs$locReg)=="down_stream"),"correlation"])),col="blue")
points(as.numeric(as.character(eQTLs[which(as.character(eQTLs$locReg)=="up_stream"),"distance"])),as.numeric(as.character(eQTLs[which(as.character(eQTLs$locReg)=="up_stream"),"correlation"])),col="green")

legend("topright",legend = c("up_stream","down_stream"),col=c("green","blue"),pch = 1)


tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[!is.na(tmp$correlation),]
misann <- tmp[which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),"locReg"]
inIR <- tmp[which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),"locReg"]
tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
novelExons <- tmp[!is.na(tmp$correlation),"locReg"]

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- tmp[!is.na(tmp$correlation),"locReg"]


rbind(misann=table(misann),inIR=table(inIR),novelExons=table(novelExons),unkwon=table(tmp),
      NoCor=table(eQTLs[is.na(eQTLs$correlation),"locReg"]),
      total=table(eQTLs[which(as.numeric(as.character(eQTLs$uniquePortion))>0),"locReg"]))

rm(tmp)



table(eQTLs[which(as.numeric(as.character(eQTLs$uniquePortion))>0),"locReg"])



#### we now look at the expression of the intergenic regions at the three different categories

load("data/results/novelIntergenicRegionsFinal/finalNovelIntergeniceQTL.PUTM.rda")

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[!is.na(tmp$correlation),]
misann <- tmp[which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),"region"]
inIR <- tmp[which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),"region"]
tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
novelExons <- tmp[!is.na(tmp$correlation),"region"]

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- tmp[!is.na(tmp$correlation),"region"]

load("data/expr/rawCounts/intergenic/exprIntergenic.PUTM.rda")
exprIntergenic <- exprIntergenic[,5:ncol(exprIntergenic)]

par(mar=c(3, 5, 5, 3))
boxplot(cbind(misAnnotation=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(misann)),]))),
              indIntergenic=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(inIR)),]))),
              novelExon=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(novelExons)),]))),
              unknown=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(tmp)),])))),col=c("red","green","blue","grey"),main="Raw counts intergenic regions",ylab="log(counts+1)")


boxplot(cbind(misAnnotation=log1p(apply(as.matrix(exprIntergenic[unique(as.character(misann)),]),1,mean)),
              indIntergenic=log1p(apply(as.matrix(exprIntergenic[unique(as.character(inIR)),]),1,mean)),
              novelExon=log1p(apply(as.matrix(exprIntergenic[unique(as.character(novelExons)),]),1,mean)),
              unknown=log1p(apply(as.matrix(exprIntergenic[unique(as.character(tmp)),]),1,mean))),col=c("red","green","blue","grey"),main="Mean by region of raw counts",ylab="log(counts+1)")




load("data/results/novelIntergenicRegionsFinal/finalNovelIntergeniceQTL.PUTM.rda")

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[!is.na(tmp$correlation),]
misann <- tmp[which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),"exon"]
inIR <- tmp[which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),"exon"]
tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
novelExons <- tmp[!is.na(tmp$correlation),"exon"]

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- tmp[!is.na(tmp$correlation),"exon"]

load("data/expr/rawCounts/genic/exons.rda")
exprExons <- countsTable 
rm(countsTable)

par(mar=c(3, 5, 5, 3))
boxplot(cbind(misAnnotation=log1p(as.numeric(unlist(exprExons[unique(as.character(misann)),]))),
              indIntergenic=log1p(as.numeric(unlist(exprExons[unique(as.character(inIR)),]))),
              novelExon=log1p(as.numeric(unlist(exprExons[unique(as.character(novelExons)),]))),
              unknown=log1p(as.numeric(unlist(exprExons[unique(as.character(tmp)),])))),col=c("red","green","blue","grey"),main="Raw counts exons",ylab="log(counts+1)")


boxplot(cbind(misAnnotation=log1p(apply(as.matrix(exprExons[unique(as.character(misann)),]),1,mean)),
              indIntergenic=log1p(apply(as.matrix(exprExons[unique(as.character(inIR)),]),1,mean)),
              novelExon=log1p(apply(as.matrix(exprExons[unique(as.character(novelExons)),]),1,mean)),
              unknown=log1p(apply(as.matrix(exprExons[unique(as.character(tmp)),]),1,mean))),col=c("red","green","blue","grey"),main="Mean by exons of raw counts",ylab="log(counts+1)")



##combined both in one plot


load("data/results/novelIntergenicRegionsFinal/finalNovelIntergeniceQTL.PUTM.rda")

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[!is.na(tmp$correlation),]
misann <- tmp[which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]
inIR <- tmp[which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- eQTLs[which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
novelExons <- tmp[!is.na(tmp$correlation),]

tmp <- eQTLs[-which(as.numeric(as.character(eQTLs$exonJunc))>0 & as.numeric(as.character(eQTLs$distance)) >0),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))<5000 & as.numeric(as.character(tmp$correlation))>0.2),]
tmp <- tmp[-which(as.numeric(as.character(tmp$distance))>5000 & as.numeric(as.character(tmp$correlation))<0.2),]
tmp <- tmp[!is.na(tmp$correlation),]

load("data/expr/rawCounts/intergenic/exprIntergenic.PUTM.rda")
exprIntergenic <- exprIntergenic[,5:ncol(exprIntergenic)]


load("data/expr/rawCounts/genic/exons.rda")
exprExons <- countsTable 
rm(countsTable)

exprAll <- as.data.frame(rbind(cbind(expr=log1p(as.numeric(unlist(exprExons[unique(as.character(misann$exon)),]))),cat="misAnn",typ="exons"),
      cbind(expr=log1p(as.numeric(unlist(exprExons[unique(as.character(inIR$exon)),]))),cat="inIR",typ="exons"),
      cbind(novelExon=log1p(as.numeric(unlist(exprExons[unique(as.character(novelExons$exon)),]))),cat="novelExon",typ="exons"),
      cbind(novelExon=log1p(as.numeric(unlist(exprExons[unique(as.character(tmp$exon)),]))),cat="unknown",typ="exons"),
      cbind(expr=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(misann$region)),]))),cat="misAnn",typ="regions"),
      cbind(expr=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(inIR$region)),]))),cat="inIR",typ="regions"),
      cbind(novelExon=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(novelExons$region)),]))),cat="novelExon",typ="regions"),
      cbind(novelExon=log1p(as.numeric(unlist(exprIntergenic[unique(as.character(tmp$region)),]))),cat="unknown",typ="regions")))

exprAll$expr <- as.numeric(as.character(exprAll$expr))




par(mar=c(8, 5, 5, 3))
boxplot(expr~typ+cat,data=exprAll, exprcol=c("red","green","blue","grey"),main="Raw counts",ylab="log(counts+1)",las=2,
        col=c("green","green","red","red","blue","blue","grey","grey"))



exprAll <- as.data.frame(rbind(cbind(expr=log1p(apply(as.matrix(exprExons[unique(as.character(misann$exon)),]),1,mean)),cat="misAnn",typ="exons"),
                               cbind(expr=log1p(apply(as.matrix(exprExons[unique(as.character(inIR$exon)),]),1,mean)),cat="inIR",typ="exons"),
                               cbind(novelExon=log1p(apply(as.matrix(exprExons[unique(as.character(novelExons$exon)),]),1,mean)),cat="novelExon",typ="exons"),
                               cbind(novelExon=log1p(apply(as.matrix(exprExons[unique(as.character(tmp$exon)),]),1,mean)),cat="unknown",typ="exons"),
                               cbind(expr=log1p(apply(as.matrix(exprIntergenic[unique(as.character(misann$region)),]),1,mean)),cat="misAnn",typ="regions"),
                               cbind(expr=log1p(apply(as.matrix(exprIntergenic[unique(as.character(inIR$region)),]),1,mean)),cat="inIR",typ="regions"),
                               cbind(novelExon=log1p(apply(as.matrix(exprIntergenic[unique(as.character(novelExons$region)),]),1,mean)),cat="novelExon",typ="regions"),
                               cbind(novelExon=log1p(apply(as.matrix(exprIntergenic[unique(as.character(tmp$region)),]),1,mean)),cat="unknown",typ="regions")))


exprAll$expr <- as.numeric(as.character(exprAll$expr))

par(mar=c(8, 5, 5, 3))
boxplot(expr~typ+cat,data=exprAll, exprcol=c("red","green","blue","grey"),main="Mean raw counts",ylab="log(counts+1)",las=2,
        col=c("green","green","red","red","blue","blue","grey","grey"))





par(mar=c(3, 5, 5, 3))
boxplot(cbind(misAnnotation=log1p(as.numeric(unlist(exprExons[unique(as.character(misann)),]))),
              indIntergenic=log1p(as.numeric(unlist(exprExons[unique(as.character(inIR)),]))),
              novelExon=log1p(as.numeric(unlist(exprExons[unique(as.character(novelExons)),]))),
              unknown=log1p(as.numeric(unlist(exprExons[unique(as.character(tmp)),])))),col=c("red","green","blue","grey"),main="Raw counts exons",ylab="log(counts+1)")


boxplot(cbind(misAnnotation=log1p(apply(as.matrix(exprExons[unique(as.character(misann)),]),1,mean)),
              indIntergenic=log1p(apply(as.matrix(exprExons[unique(as.character(inIR)),]),1,mean)),
              novelExon=log1p(apply(as.matrix(exprExons[unique(as.character(novelExons)),]),1,mean)),
              unknown=log1p(apply(as.matrix(exprExons[unique(as.character(tmp)),]),1,mean))),col=c("red","green","blue","grey"),main="Mean by exons of raw counts",ylab="log(counts+1)")






