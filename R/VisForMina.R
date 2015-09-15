

library("biomaRt")
## we select the archive version
listMarts(host="Jun2013.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

load("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr17.rda")
load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]


chr <- 17
start <- 44574341   
end <- 44575417
gen = "hg19"


fullCovtmp <- fullCov$chr17[(start-1000):(end+10000),PUTM$A.CEL_file]
dim(fullCovtmp)
rownames(fullCovtmp) <- (start-1000):(end+10000)
     
tail(data.frame(fullCovtmp))
     
filters <- c("chromosomal_region")


chromoReg <- paste0(chr,":",start-1000,":",end+10000,":-1,",
                    chr,":",start-1000,":",end+10000,":1")

     
annReg <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                                "exon_chrom_end","strand","gene_biotype",
                                "ensembl_gene_id",'ensembl_exon_id',
                                'ensembl_transcript_id',"external_gene_id")
                   , filters=filters, values=list(chromoReg), mart=ensembl)


width <- (annReg$exon_chrom_end - annReg$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
annReg <- cbind(annReg[,1:3],width,annReg[,4:9])  
rm(width)

annReg$strand <- gsub("-1","-",annReg$strand)
annReg$strand <- gsub("1","+",annReg$strand)
annReg$strand <- as.factor(annReg$strand)
## we now change the column names

colnames(annReg) <- c("chromosome","start","end","width","strand",
                         "feature","gene","exon","transcript","symbol")  



grtrack <- GeneRegionTrack(annReg, genome = gen,
                           chromosome = chr)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

## we load the dat
meanCov <- apply(data.frame(fullCovtmp),1,mean)
dat <- data.frame('chr17',as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome="chr17",genome="hg19",type="histogram")
dtrack_heat <- DataTrack(range=data_g,chromosome="chr17",genome="hg19",type="heatmap")

# plotTracks(list(gtrack, dtrack,grtrack),transcriptAnnotation = "symbol")
# 
plotTracks(list(gtrack,dtrack,grtrack,dtrack_heat), from = (start-1000),
            to = (end+10000), chromosome = chr, 
            type = "histogram",transcriptAnnotation = "symbol")


plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "symbol")






library("biomaRt")
## we select the archive version
listMarts(host="Jun2013.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

load("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr9.rda")
load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]


chr <- 9
start <- 3146055   
end <- 3146576
gen = "hg19"


fullCovtmp <- fullCov$chr9[(start-70000):(end+70000),PUTM$A.CEL_file]
dim(fullCovtmp)
rownames(fullCovtmp) <- (start-70000):(end+70000)

tail(data.frame(fullCovtmp))

filters <- c("chromosomal_region")


chromoReg <- paste0(chr,":",start-70000,":",end+70000,":-1,",
                    chr,":",start-70000,":",end+70000,":1")


annReg <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=list(chromoReg), mart=ensembl)


width <- (annReg$exon_chrom_end - annReg$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
annReg <- cbind(annReg[,1:3],width,annReg[,4:9])  
rm(width)

annReg$strand <- gsub("-1","-",annReg$strand)
annReg$strand <- gsub("1","+",annReg$strand)
annReg$strand <- as.factor(annReg$strand)
## we now change the column names

colnames(annReg) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  



grtrack <- GeneRegionTrack(annReg, genome = gen,
                           chromosome = chr)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

## we load the dat
meanCov <- apply(data.frame(fullCovtmp),1,mean)
dat <- data.frame('chr9',as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome="chr9",genome="hg19",type="histogram")
dtrack_heat <- DataTrack(range=data_g,chromosome="chr9",genome="hg19",type="heatmap")

# plotTracks(list(gtrack, dtrack,grtrack),transcriptAnnotation = "symbol")
# 
plotTracks(list(dtrack,grtrack), from = (start-70000),
           to = (end+70000), chromosome = chr, 
           type = "histogram",transcriptAnnotation = "symbol")


plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "symbol")




## we select the archive version
listMarts(host="Jun2013.archive.ensembl.org")
## we select the database from archive 72
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


load("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr1.rda")
load("data/general/sampleInfo.rda")
SNIG <- sampleInfo[which(sampleInfo$U.Region_simplified=="SNIG"),]

chr <- 1

start <- 232638989   
end <- 232739836
gen = "hg19"


fullCovtmp <- fullCov$chr1[(start-1000):(end+10000),SNIG$A.CEL_file]
dim(fullCovtmp)
rownames(fullCovtmp) <- (start-1000):(end+10000)

tail(data.frame(fullCovtmp))

filters <- c("chromosomal_region")


chromoReg <- paste0(chr,":",start-1000,":",end+10000,":-1,",
                    chr,":",start-1000,":",end+10000,":1")


annReg <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                             "exon_chrom_end","strand","gene_biotype",
                             "ensembl_gene_id",'ensembl_exon_id',
                             'ensembl_transcript_id',"external_gene_id")
                , filters=filters, values=list(chromoReg), mart=ensembl)


width <- (annReg$exon_chrom_end - annReg$exon_chrom_start) +1

## the strand field need to be converted
## convert 1 to + and -1 to -
annReg <- cbind(annReg[,1:3],width,annReg[,4:9])  
rm(width)

annReg$strand <- gsub("-1","-",annReg$strand)
annReg$strand <- gsub("1","+",annReg$strand)
annReg$strand <- as.factor(annReg$strand)
## we now change the column names

colnames(annReg) <- c("chromosome","start","end","width","strand",
                      "feature","gene","exon","transcript","symbol")  

library(Gviz)

grtrack <- GeneRegionTrack(annReg, genome = gen,
                           chromosome = chr)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

##load(file="/home/seb/")

## we load the dat
meanCov <- apply(data.frame(fullCovtmp),1,mean)
dat <- data.frame('chr1',as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
rm(meanCov)
colnames(dat) <- c("chr","start","end","coverage")
data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
dtrack <- DataTrack(range=data_g,chromosome="chr1",genome="hg19",type="histogram")
dtrack_heat <- DataTrack(range=data_g,chromosome="chr1",genome="hg19",type="heatmap")

# plotTracks(list(gtrack, dtrack,grtrack),transcriptAnnotation = "symbol")
# 
plotTracks(list(dtrack,grtrack), from = (start-1000),
           to = (end+10000), chromosome = chr, 
           type = "histogram",transcriptAnnotation = "symbol")


plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "symbol")



65505809 65506637





