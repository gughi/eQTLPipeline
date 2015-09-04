


setwd("/home/guelfi/eQTLPipeline")
sink("logExonicIntornic.log")
nCores <- 15
cat(paste("Number of cores",nCores,"\n"))
library(devtools)
load_all()

cat(paste("Processing PUTM region \n"))

## Intergenic data was counts were generated using DerFinder
load("data/expr/rawCounts/intergenic/PUTM/allChromosomes.rda")
## load data has 2 different objects: annotation,coverage
## load the sample info
load("data/general/sampleInfo.rda")
## we select only samples from PUTM

PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]
rm(sampleInfo)

exprIntergenic <- cbind(annotation$seqnames,annotation$start,annotation$end,annotation$width,coverage)
colnames(exprIntergenic) <- c("chr","start","end","width",colnames(coverage))
rm(annotation,coverage)    


## select regions with length => 100bp
exprIntergenic <- exprIntergenic[which(exprIntergenic$width >= 100),]
# We don't do filtering since it doesn't make any since for DERFINDER; 
# regions are detected based on the expression
# creation of identifiers for the regions
IDs <- paste0("DER",c(1:nrow(exprIntergenic)))
rownames(exprIntergenic) <- IDs
rm(IDs)

cat(paste("number of intergenic regions >100bp",nrow(exprIntergenic)))



# load the library size
librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

# convert in RPKM
library(easyRNASeq)

# load the GC content intergenic and region length        
length <- exprIntergenic$width    
names(length) <- rownames(exprIntergenic)
stopifnot(identical(colnames(as.matrix(exprIntergenic[,as.character(PUTM$A.CEL_file)]))
                    ,names(librarySize)))
stopifnot(identical(rownames(exprIntergenic),names(length)))              

#convert in RPKM
RPKM.std <- RPKM(as.matrix(exprIntergenic[,as.character(PUTM$A.CEL_file)])
                 , NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
genesList <- rownames(RPKM.std) 
rm(RPKM.std,length)
exprIntergenic <- exprIntergenic[as.character(genesList),]

cat(paste("number of regions included in the analysis:",nrow(exprIntergenic)))

cat("Calculating the GC content")
GCcontent <- GCcalculation(exprIntergenic[,1:4],genRef="/home/ukbec/bowtie2Index/genome37.72.fa"
                           ,pathBedtools = "/apps/BEDTools/2.24.0/bin/bedtools")
GCcontent <- as.data.frame(cbind(GCcontent[as.character(rownames(exprIntergenic)),],
                                 exprIntergenic$width))    
rownames(GCcontent) <- rownames(exprIntergenic)
colnames(GCcontent) <- c("GCcontent","length")



