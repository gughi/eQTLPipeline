

load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.SNIG")

exonAnn <- read.delim("data/general/exongrps.log",sep=" ",skip=5,header=F)
head(exonExonJun[,1:5])
head(exonAnn)


library("biomaRt")
## we select the archive version
listMarts(host="Jun2013.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

filters <- c("ensembl_gene_id")

geneID <- "ENSG00000188906" 
gen = "hg19"

geneAnn <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=geneID, mart=ensembl)


width <- (geneAnn$exon_chrom_end - geneAnn$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
geneAnn <- cbind(geneAnn[,1:3],width,geneAnn[,4:9])  
rm(width)

geneAnn$strand <- gsub("-1","-",geneAnn$strand)
geneAnn$strand <- gsub("1","+",geneAnn$strand)
geneAnn$strand <- as.factor(geneAnn$strand)
## we now change the column names

colnames(geneAnn) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  


## select the interval to get the exon-exon junctions
exonExonGene <- exonExonJun[which(exonExonJun$Chr == unique(geneAnn$chromosome) 
                  & exonExonJun$TSS >= (min(geneAnn$start))
                  & exonExonJun$TSS <= (max(geneAnn$end))),]




exonExonGene <- RPKM.cqn[which(exonExonJun$Chr == unique(geneAnn$chromosome) 
                                  & exonExonJun$TSS >= (min(geneAnn$start))
                                  & exonExonJun$TSS <= (max(geneAnn$end))),]




exonAnnbyGene <- exonAnn[grep(geneID,exonAnn$V6),]

head(exonAnnbyGene[,1:5])
head(exonExonGene[,1:4])



head(exonAnnbyGene)

library(Gviz)

grtrack <- GeneRegionTrack(geneAnn, genome = gen,
                           chromosome = unique(geneAnn$chromosome), name = geneID)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = unique(geneAnn$chromosome))

##load(file="/home/seb/")

## we load the dat

meanCov <- apply(data.frame(exonExonGene[,-c(1:4)]),1,mean)

head(meanCov)
dat <- data.frame(chr=NA, start=NA,end=NA,coverage=NA)[numeric(0), ]
for(i in 2:nrow(exonExonGene))
  {
  
      dat <- rbind(dat,test(exonExonGene[i,1],i))
  }

for(i in 2:nrow(exonExonGene))
{
  
  dat <- rbind(dat,test(exonExonGene[i,2],i))
}


  
test <- function(x,i)
{
  tmp <- exonAnnbyGene[which(exonAnnbyGene$V5 == x),]
  tmp <- data.frame(paste0("chr",tmp[,2]),
                tmp[,3],tmp[,4],meanCov[i])
  colnames(tmp) <- c("chr","start","end","coverage")
  return(tmp)
}

## dat <- data.frame(paste0("chr",unique(geneAnn$chromosome)),as.numeric(exonExonGene$TSS),as.numeric(exonExonGene$TSS),as.vector(meanCov))                     
rm(meanCov)
class(dat$start)


data_g <- with(dat, GRanges(chr, IRanges(as.numeric(as.character(start)),as.numeric(as.character(end))), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome=paste0("chr",unique(geneAnn$chromosome)),genome="hg19",name=geneID,type="histogram")
dtrack_heat <- DataTrack(range=data_g,chromosome=paste0("chr",unique(geneAnn$chromosome)),genome="hg19",name=geneID,type="heatmap")

## itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "transcript", type=c("coverage","sashimi"))

dtrack





traID <- "ENST00000298910"



exonAnnbyTra <- exonAnn[grep(geneID,exonAnn$V6),]


meanCov <- apply(data.frame(exonExonGene[,-c(1:4)]),1,mean)

head(meanCov)
dat <- data.frame(chr=NA, start=NA,end=NA,coverage=NA)[numeric(0), ]
for(i in 2:nrow(exonExonGene))
{
  
  dat <- rbind(dat,test(exonExonGene[i,1],i))
}

for(i in 2:nrow(exonExonGene))
{
  
  dat <- rbind(dat,test(exonExonGene[i,2],i))
}



test <- function(x,i)
{
  tmp <- exonAnnbyGene[which(exonAnnbyGene$V5 == x),]
  tmp <- data.frame(paste0("chr",tmp[,2]),
                    tmp[,3],tmp[,4],meanCov[i])
  colnames(tmp) <- c("chr","start","end","coverage")
  return(tmp)
}



## try to get the sashimi plots

# source("https://bioconductor.org/biocLite.R")
# biocLite("ggbio")



# devtools::install_github("pkimes/spliceclust")
# 
# install.packages("Rcpp")
# library(Rcpp)



library(spliceclust)

load("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr12.rda")


tt <- apply(fullCov$[exonAnnbyGene$V3[2]:exonAnnbyGene$V4[2],],1,mean)



## For Daniah
geneID <- "ENSG00000188906" 

load("data/expr/normalisedCounts/genic/exonExonJunc/RPKM.cqn.SNIG")

exonAnn <- read.delim("data/general/exongrps.log",sep=" ",skip=5,header=F)
head(exonAnn)


LRRK2ExExJunc <- exonAnn[grep(geneID,exonAnn$V6),]

head(LRRK2ExExJunc)

tmp <- expand.grid(LRRK2ExExJunc$V5,LRRK2ExExJunc$V5)
tmp <- paste(tmp$Var1,tmp$Var2,sep="_")
tmp <- intersect(rownames(RPKM.cqn),tmp)


LRRK2Junc <- RPKM.cqn[tmp,]


write.csv(LRRK2Junc,file="C:/Users/mguelfi/Desktop/Projects/daniah/exonExonJunc.csv")
head(exonAnn)
LRRK2ExExJunc <- LRRK2ExExJunc[,c(2:5)]
colnames(LRRK2ExExJunc) <- c("chr","start","end","exonID")
write.csv(LRRK2ExExJunc,file="C:/Users/mguelfi/Desktop/Projects/daniah/AnnLRRK2.csv")


library(GenomicFeatures)
> gff_file <- system.file("extdata", "GFF3_files", "a.gff3",
                          +                         package="GenomicFeatures")
txdb <- makeTranscriptDbFromGFF("/home/seb/forDaniah/LRRK2.gtf", format="gtf")

tmp <- as.data.frame(exonsBy(txdb, by="gene"))
tmp <- paste(tmp$seqnames,tmp$start,tmp$end,tmp$group_name,sep="\t")


write.table(data.frame(tmp), file = "/home/seb/forDaniah/LRRK2.BED", row.names = F, 
            col.names = F, quote = F)
rm(tmp)

cmd <- paste0("bedtools getfasta -fi /home/seb/reference/genome37.72.fa  -bed /home/seb/forDaniah/LRRK2.BED -fo /home/seb/forDaniah/sequences.out")

head(read.delim("/home/seb/forDaniah/sequences.out",header=F))










