


exonExonJun <- read.delim("data/expr/rawCounts/genic/allSamples.2norm.merged")
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















