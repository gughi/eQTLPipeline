## Checking the questions Jana asked about MAPT


MAPT <- "ENSG00000186868"

## Is it in both PUTM and SNIG
library("biomaRt")

## this is what we are seeing for the sentinal SNPs

eQTLPUTM <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

genes <- unlist(lapply(strsplit(as.character(eQTLPUTM$gene),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))


eQTLPUTM[which(genes %in% MAPT),]


#             snps                 gene     tstat       pvalue          FDR       beta        myFDR degree
# 1809       chr17:43757161 ENSG00000186868:E023 -7.408543 4.149581e-11 4.808721e-10 -0.5058429 2.342024e-07      1
# 6306  chr17:44212211:C_CA ENSG00000186868:E012 -6.320799 7.310657e-09 4.938189e-08 -0.5424641 4.126135e-05      1
# 7115       chr17:44359663 ENSG00000186868:E016 -7.437221 3.608864e-11 1.083882e-09 -0.4931524 2.036843e-07      1
# 7948       chr17:44304884 ENSG00000186868:E022 -6.887455 5.107501e-10 5.320318e-09 -0.3084519 2.882674e-06      1
# 9560       chr17:44308547 ENSG00000186868:E007 -7.000464 2.976946e-10 9.073704e-09 -0.3777567 1.680189e-06      1
# 9982       chr17:43692338 ENSG00000186868:E015 -6.564069 2.355521e-09 3.537722e-08 -0.3037446 1.329456e-05      1
# 10384      chr17:44359663 ENSG00000186868:E017 -7.397417 4.380386e-11 1.474154e-09 -0.4646552 2.472290e-07      1
# 11978      chr17:44351600 ENSG00000186868:E009 -7.207908 1.097629e-10 7.403821e-09 -0.4020811 6.195017e-07      1



eQTLSNIG <- read.delim("data/results/finaleQTLs/exons.SNIG.txt",sep=" ")

genes <- unlist(lapply(strsplit(as.character(eQTLSNIG$gene),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))


eQTLSNIG[which(genes %in% MAPT),]

## rs183211 chr17:44788310


rm(eQTLSNIG,genes,eQTLPUTM)
## Unsentinal SNPs

eQTLPUTM   <- read.delim(file="data/results/finaleQTLs/exons.unsentinalised.PUTM.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLPUTM$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
head(df)
eQTLPUTM <- df
save(eQTLPUTM,file= "data/results/finaleQTLs/exons.unsentinalised.formatted.PUTM.txt")

genes <- unlist(lapply(strsplit(as.character(eQTLPUTM[,2]),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

write.csv(unique(eQTLPUTM[which(genes %in% MAPT),1]),file="tmp/SNPs.MAPT.PUTM.forJANA.csv")

rm(genes,df,tmp,eQTLPUTM)

eQTLSNIG   <- read.delim(file="data/results/finaleQTLs/exons.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLSNIG$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
##head(df)
rownames(df) <- NULL
head(df)
eQTLSNIG <- df
save(eQTLSNIG,file= "data/results/finaleQTLs/exons.unsentinalised.formatted.SNIG.txt")

genes <- unlist(lapply(strsplit(as.character(eQTLSNIG[,2]),":"),function(x){x[1]}))
genes <- unlist(sapply(strsplit(as.character(genes),"+",fixed=T),function(x){x[1]}))

## there are not exons regulated in MAPT
##write.csv(unique(eQTLSNIG[which(genes %in% MAPT),1]),file="tmp/SNPs.MAPT.SNIG.forJANA.csv")

rm(genes,df,tmp,eQTLSNIG)


## JV: Any exon level eQTL for exon 10?
## SG: below, the coordinates of the exons that are regulated.

exonsMAPT <- c("ENSG00000186868:E023",
  "ENSG00000186868:E012",
  "ENSG00000186868:E016",
  "ENSG00000186868:E022",
  "ENSG00000186868:E007",
  "ENSG00000186868:E015",
  "ENSG00000186868:E017",
  "ENSG00000186868:E009")


transcriptomeInfo <- read.csv("data/general/transcriptomeInfo.csv",row.names=4)
head(transcriptomeInfo)
transcriptomeInfo[exonsMAPT,]

# ENSG00000186868:E023 44067244 44067441   198
# ENSG00000186868:E012 44051751 44051837    87
# ENSG00000186868:E016 44055747 44055791    45
# ENSG00000186868:E022 44064406 44064461    56
# ENSG00000186868:E007 44039704 44039836   133
# ENSG00000186868:E015 44055741 44055746     6
# ENSG00000186868:E017 44055792 44055806    15
# ENSG00000186868:E009 44049225 44049311    87





rm(exonsMAPT,transcriptomeInfo)
MAPT <- "ENSG00000186868"

load("data/results/finaleQTLs/eQTL.ExExJun.PUTM.rda")

eQTL.PUTM[which(eQTL.PUTM$geneID %in% MAPT),]


# snps        exExID    tstat       pvalue          FDR      beta        myFDR degree          geneID
# 236519_236521       chr17:44848314 236519_236521 4.851766 4.508251e-06 3.608964e-04 0.3926621 2.544457e-02      1 ENSG00000186868
# 236511_236518 chr17:44782224:A_AGT 236511_236518 7.789597 6.423347e-12 8.512700e-11 0.7976702 3.625337e-08      1 ENSG00000186868
# transID chrExon1 startExon1 endExon1 chrExon2 startExon2 endExon2 differentGene distanceExons numOfTrans hasUniExo
# ENST00000415613;ENST00000571311       17   44039703 44039836       17   44049224 44049311         FALSE          9388          1     FALSE
# ENST00000570299;ENST00000571311       17   43971941 43972052       17   44039686 44039836         FALSE         67634          1     FALSE


load("data/results/finaleQTLs/eQTL.ExExJun.SNIG.rda")
eQTL.SNIG[which(eQTL.SNIG$geneID %in% MAPT),]

# snps        exExID    tstat       pvalue          FDR     beta       myFDR degree          geneID                         transID
# 236511_236518 chr17:44034575 236511_236518 5.460145 9.587868e-07 3.357296e-05 1.143994 0.005411393      1 ENSG00000186868 ENST00000570299;ENST00000571311
# chrExon1 startExon1 endExon1 chrExon2 startExon2 endExon2 differentGene distanceExons numOfTrans hasUniExo external_gene_id   gene_biotype
# 236511_236518       17   43971941 43972052       17   44039686 44039836         FALSE         67634          1     FALSE             MAPT protein_coding


load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTM[which(eQTLPUTM$gene %in% MAPT),]


load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
eQTLSNIG[which(eQTLSNIG$gene %in% MAPT),]

