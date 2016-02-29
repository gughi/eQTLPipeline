
## load 
setwd("C:/Users/mguelfi/projectsR/spliceclust-master/")
library(devtools)
load_all()
##
setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()

## load the exon-exon junctions
library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

## we load counts and mapping file for the exon counts
load("data/expr/rawCounts/genic/fullExExJun.rda")

MAPT <- "ENSG00000186868"

## we select the exon junctions for MAPT
MAPTexonJunc <- mapExon[grep(MAPT,mapExon$geneID),]

idx <- unique(c(which(expr$Exon1ID %in% MAPTexonJunc$exonID),which(expr$Exon2ID %in% MAPTexonJunc$exonID)))

exprGene <- expr[idx,]


library(GenomicRanges)

head(MAPTexonJunc)
rm(idx)


colnames(exprGene) <- gsub("Sample_A653_","s",colnames(exprGene))

head(exprGene)

juncPosition <- NULL
for(i in 1:nrow(exprGene))
{
  if (is.element(exprGene$Exon1ID[i],MAPTexonJunc$exonID) & is.element(exprGene$Exon2ID[i],MAPTexonJunc$exonID))
  {
  juncPosition <- rbind(juncPosition,c(chr=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon1ID[i]),"chr"]),  
       start=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon1ID[i]),"end"]),
       end=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon2ID[i]),"start"]),kind="j",
       ##data.frame(t(setNames(rep(0,10),colnames(exprGene)[5:14])))
       exprGene[i,5:ncol(exprGene)]
       ))
  # we add also the exons that are needed for the plot
  juncPosition <- rbind(juncPosition,c(chr=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon1ID[i]),"chr"]),  
                                       start=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon1ID[i]),"start"]),
                                       end=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon1ID[i]),"end"]),kind="e",
                                       #exprGene[i,5:24]
                                       data.frame(t(setNames(rep(0,ncol(exprGene)-4),colnames(exprGene)[5:ncol(exprGene)])))
                                       ))
  
  juncPosition <- rbind(juncPosition,c(chr=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon1ID[i]),"chr"]),  
                                       start=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon2ID[i]),"start"]),
                                       end=as.character(MAPTexonJunc[which(MAPTexonJunc$exonID %in% exprGene$Exon2ID[i]),"end"]),kind="e",
                                       #exprGene[i,5:24]
                                       data.frame(t(setNames(rep(0,ncol(exprGene)-4),colnames(exprGene)[5:ncol(exprGene)])))
                                       ))
  }
}



juncPosition <- apply(juncPosition,2,unlist)
juncPosition <- data.frame(unique(juncPosition),stringsAsFactors =  F )
juncPosition[,1:3] <- apply(juncPosition[,1:3], 2,as.numeric)
juncPosition[,5:(ncol(exprGene))] <- apply(juncPosition[,5:(ncol(exprGene))], 2,as.numeric)

# juncPosition <- juncPosition[1:3,]
# juncPosition$names <- 1:3
# juncPosition$gIdx <- "test"
# juncPosition$gStart <- 43972052
# juncPosition$gStop <-  44039836 


gr <- makeGRangesFromDataFrame(juncPosition,
                               keep.extra.columns = T )


grl <- split(gr, mcols(gr)$kind)
gr_cc <- concomp(grl)

splicegrahm(gr_cc)
splicegrahm(gr_cc,j_incl = TRUE)



library("biomaRt")
## we select the archive version
listMarts(host="Jun2013.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


filters <- c("ensembl_gene_id")


annReg <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand",
                             'ensembl_exon_id',
                             'ensembl_transcript_id',"rank")
                , filters=filters, values=list(MAPT), mart=ensembl)


width <- (annReg$exon_chrom_end - annReg$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
annReg <- cbind(annReg[,1:3],width,annReg[,4:7])  
rm(width)

annReg$chromosome <- paste0("chr",annReg$chromosome)
annReg$strand <- gsub("-1","-",annReg$strand)
annReg$strand <- gsub("1","+",annReg$strand)
annReg$ensembl_exon_id <- as.numeric(gsub("ENSE","",annReg$ensembl_exon_id))
annReg$ensembl_transcript_id <- as.numeric(gsub("ENST","",annReg$ensembl_transcript_id))
annReg$strand <- as.factor(annReg$strand)
## we now change the column names

head(annReg)

colnames(annReg) <- c("chromosome","start","end","width","strand",
                      "exon_id","transcript","exon_rank")


ann <- makeGRangesFromDataFrame(annReg,
                               keep.extra.columns = T )

ann <-split(ann, mcols(ann)$transcript)


splicegrahm(gr_cc, txlist = ann)






suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)[paste0("chr", 1:22)] <- TRUE
exbytx <- exonsBy(txdb, "tx")



## working

tmp <- data.frame(rbind(c(51532348,51532468,121,33127),c(51532468,51532598,131,33128),
      c(51532469,51532597,129,33129),c(51532598,51532713,116,33130),
      c(51532713, 51534044,  1332, 33131)),stringsAsFactors = F)

colnames(tmp) <-c("start","end","width","names")
tmp$kind <- c("e", "j", "e", "e", "j")

tmp$chr <- 9
tmp$gIdx <- "gene9317"         
tmp$gStart <- 51532348
tmp$gStop <- 51538261 
tmp$s1 <- c(2.60331,0.00000,2.59690,2.15517,0.00000)
tmp$s2  <- c(5.32231,0.00000,10.37210,7.28448,0.00000)        



gr1 <- makeGRangesFromDataFrame(tmp,
                               keep.extra.columns = T )

mcols(gr)[1:3, 1:6]


gl1 <- split(gr1, mcols(gr1)$kind)

cc1 <- concomp(gl1)
splicegrahm(cc1)

splicegrahm(cc1,txlist = exbytx)















