### Manuel Sebastian Guelfi
### 20/4/15
### Script to visualise transcripts 
 
## source("http://bioconductor.org/biocLite.R")
## biocLite("Gviz")
## source("http://bioconductor.org/biocLite.R")
## biocLite("biomaRt")
## needed also to update this library to make itracker work
##biocLite("rtracklayer")


library("biomaRt")
## we select the archive version
listMarts(host="Jun2013.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")
## we select the the homo_sapiens database
##ensembl <-  useDataset("hsapiens_gene_ensembl",mart=ensembl)

##filters <- listFilters(ensembl)
##attributes <- listAttributes(ensembl)

## we select the lincRNA
filters <- c("chromosome_name","biotype","status")
chr <- "21"
geneType <- "lincRNA"
sta <- "KNOWN"
lincRNA21 <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                                "exon_chrom_end","strand","gene_biotype",
                                "ensembl_gene_id",'ensembl_exon_id',
                                'ensembl_transcript_id',"external_gene_id")
              , filters=filters, values=list(chr,geneType,sta), mart=ensembl)

rm(geneType,sta,filters,chr)
width <- (lincRNA21$exon_chrom_end - lincRNA21$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
lincRNA21 <- cbind(lincRNA21[,1:3],width,lincRNA21[,4:9])  
rm(width)

lincRNA21$strand <- gsub("-1","-",lincRNA21$strand)
lincRNA21$strand <- gsub("1","+",lincRNA21$strand)
lincRNA21$strand <- as.factor(lincRNA21$strand)
## we now change the column names

colnames(lincRNA21) <- c("chromosome","start","end","width","strand",
  "feature","gene","exon","transcript","symbol")  

## load the library
library(Gviz)

load("/home/seb/DERFINDER/fullCoverageChr21.rda")

## we load the data for the example
chr = 21
gen = "hg19"
gene <- "ENSG00000215386"
lincRNA21tmp <- lincRNA21[ which(lincRNA21$gene %in% gene),]


### Select of the bp for the gene considered

startStop <- getBM(attributes=c("start_position","end_position"), filters="ensembl_gene_id", values=list(gene), mart=ensembl)

fullCovtmp <- fullCov$chr21[(startStop$start_position-5000):(startStop$end_position+5000),]

rownames(fullCovtmp) <- (startStop$start_position-5000):(startStop$end_position+5000)

tail(data.frame(fullCovtmp))



grtrack <- GeneRegionTrack(lincRNA21tmp, genome = gen,
                          chromosome = chr, name = gene)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
      
##load(file="/home/seb/")

## we load the dat
meanCov <- apply(data.frame(fullCovtmp),1,mean)
dat <- data.frame('chr21',as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome="chr21",genome="hg19",name=gene,type="histogram")
dtrack_heat <- DataTrack(range=data_g,chromosome="chr21",genome="hg19",name=gene,type="heatmap")

## itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "transcript")                     

### zoom in

plotTracks(c(dtrack,grtrack,dtrack_heat ), from = startStop$start_position-5000,
           to = 17450000, chromosome = "chr21")



## to get the annotation from ENSEMBL
##biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
##                                    chromosome = chr, 
##                                    start = 20000000, end = 21000000,
##                                    name = "ENSEMBL")
##plotTracks(biomTrack)


############################################################
## Novel LincRNA found in the version 75: ENSG00000273492 ##
############################################################



library("biomaRt")
## we select the archive version
listMarts(host="feb2014.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="feb2014.archive.ensembl.org")
## we select the the homo_sapiens database
ensembl <-  useDataset("hsapiens_gene_ensembl",mart=ensembl)



##filters <- listFilters(ensembl)
##attributes <- listAttributes(ensembl)

## we select the lincRNA
filters <- c("chromosomal_region")


gen = "hg19"
chr <- "21"
geneType <- "lincRNA"
gene <- "ENSG00000273492"

startStop <- getBM(attributes=c("start_position","end_position"), filters="ensembl_gene_id", values=list(gene), mart=ensembl)

chromoReg <- paste0(chr,":",startStop$start_position,":",startStop$end_position,":-1,",
                    chr,":",startStop$start_position,":",startStop$end_position,":1")

lincRNA21 <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                                "exon_chrom_end","strand","gene_biotype",
                                "ensembl_gene_id",'ensembl_exon_id',
                                'ensembl_transcript_id',"external_gene_id")
                   , filters=filters, values=list(chromoReg), mart=ensembl)

rm(geneType,sta,filters)
width <- (lincRNA21$exon_chrom_end - lincRNA21$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
lincRNA21 <- cbind(lincRNA21[,1:3],width,lincRNA21[,4:9])  
rm(width)

lincRNA21$strand <- gsub("-1","-",lincRNA21$strand)
lincRNA21$strand <- gsub("1","+",lincRNA21$strand)
lincRNA21$strand <- as.factor(lincRNA21$strand)
## we now change the column names

colnames(lincRNA21) <- c("chromosome","start","end","width","strand",
                         "feature","gene","exon","transcript","symbol")  

## load the library
library(Gviz)

load("/home/seb/DERFINDER/fullCoverageChr21.rda")

## we load the data for the example

### Select of the bp for the gene considered

fullCovtmp <- fullCov$chr21[(startStop$start_position-10000):(startStop$end_position+10000),]

rownames(fullCovtmp) <- (startStop$start_position-10000):(startStop$end_position+10000)

tail(data.frame(fullCovtmp))



grtrack <- GeneRegionTrack(lincRNA21, genome = gen,
                           chromosome = chr, name = gene)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

## we load the dat
meanCov <- apply(data.frame(fullCovtmp),1,mean)
dat <- data.frame('chr21',as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome="chr21",genome="hg19",name=gene,type="histogram")
dtrack_heat <- DataTrack(range=data_g,chromosome="chr21",genome="hg19",name=gene,type="heatmap")

## itrack <- IdeogramTrack(genome = gen, chromosome = chr)

plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "gene")                     

zoom

plotTracks(c(dtrack,grtrack ), from = (startStop$start_position-10000),
           to = (startStop$end_position+10000), chromosome = chr, 
           type = "histogram",transcriptAnnotation = "gene")

plotTracks(c(dtrack,grtrack ), from = 27580000,
           to = 27600000, chromosome = chr, 
           type = "histogram",transcriptAnnotation = "transcript")

                     
