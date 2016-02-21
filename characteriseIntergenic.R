
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

save(eQTLs,file="data/results/finaleQTLs/betaInterIntergenicNearestGene.rda")



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


    
    
    