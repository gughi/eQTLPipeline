

## load teh eQTL table
load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR <=0.05),1:7]

## exonJunc tells in how manysamples this regions has been seeing as an exon-junction
eQTLPUTM$exonJunc <- 0
## exonJuncCounts tells how many exon junctions has been seeying in total in all the samples
eQTLPUTM$exonJuncCounts <- 0


samples <- c("A653_001_CEL", "A653_012_CEL", "A653_031_CEL", "A653_048_CEL", "A653_056_CEL", "A653_057_CEL",
             "A653_063_CEL","A653_066_CEL","A653_071_CEL","A653_098_CEL","A653_106_CEL","A653_107_CEL","A653_108_CEL",
             "A653_122_CEL")

for(i in 1:length(samples))
{
  eQTLPUTM <- getSummaryNovExon(eQTLPUTM,intNovelExonPath=paste0("data/junctions/PUTM/",samples[i],"_intergenicNovelExon.rda"))
}

table(eQTLPUTM$exonJunc)

eQTLPUTM[which(eQTLPUTM$exonJunc>1),]


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


load("data/general/defIntergenic.rda")

ExonExtention <- regInformation[which(as.numeric(as.character(regInformation$overlap)) == 1),]

load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR <=0.05),1:7]

regInformation <- regInformation[as.character(eQTLPUTM$gene),]

regInformation <- regInformation[which(as.numeric(as.character(regInformation$uniquePortion))>0), ]


nrow(regInformation[which(as.numeric(as.character(regInformation$distance))<0), ])


load("data/expr/normalisedCounts/genic/geneExons/resids.PUTM.rda")
exprExonic <- resids
load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
exprIntronic <- resids
rm(resids,listNonOve)

eQTLs <- eQTLs[which(as.character(eQTLs$gene) %in% as.character(colnames(exprIntronic))),]
eQTLs <- eQTLs[which(as.character(eQTLs$gene) %in% as.character(colnames(exprExonic))),]

## remove eQTLs that are duplicated
eQTLs <- eQTLs[-which(duplicated(eQTLs[,1:2])),]

res <- NULL
for (j in 1:nrow(eQTLs))
{
  ## load the expression for each individual gene for exonic and intronic
  exprE <- exprExonic[,as.character(eQTLs$gene[j])]
  exprI <- exprIntronic[,as.character(eQTLs$gene[j])]
  
  ## add to the expression the tissue information or the exonic/intronic information
  expr <- cbind(exprE,"Exonic")
  exprtmp <- cbind(exprI,"Intronic")
  expr <- rbind(expr,exprtmp)
  rm(exprtmp,exprE,exprI)
  
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
  
  f <- betaInteraction(expr,snp0)
  ## collect summary statisctics from gene exonic
  cmd <- paste0("grep -w '",as.character(eQTLs$snps[j]),"' data/results/genic/geneExons/fullResults/PUTM/",as.character(eQTLs$gene[j]))
  statExonic <- read.delim(pipe(cmd),header=F)
  rm(cmd)
  colnames(statExonic) <- c("ge.SNP","ge.gene","ge.beta","ge.t-stat","ge.p-value","ge.FDR") 
  ## collect summary statisctics from gene intronic
  cmd <- paste0("grep -w '",as.character(eQTLs$snps[j]),"' data/results/genic/geneIntronic/fullResults/PUTM/",as.character(eQTLs$gene[j]))
  statIntronic <- read.delim(pipe(cmd),header=F)
  rm(cmd)
  colnames(statIntronic) <- c("gi.SNP","gi.gene","gi.beta","gi.t-stat","gi.p-value","gi.FDR") 
  
  restmp <- cbind(f,statExonic,statIntronic)
  colnames(restmp)[1] <- "p.interaction"
  
  if (is.null(res))
  {
    res <- restmp
  } else {
    res <- rbind(res,restmp)
  }
  rm(statExonic,statIntronic,expr,markers,markers.info,snp0,selSamples,restmp,f)
}


  




    
    
    