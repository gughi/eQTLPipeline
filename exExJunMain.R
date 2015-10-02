
setwd("/home/guelfi/eQTLPipeline")
sink("logExonExonJunctions.log")
nCores <- 15
cat(paste("Number of cores",nCores,"\n"))
library(devtools)
load_all()

## Now we correct for PEER using simple quantification Exons+Introns



expr <- read.delim("data/expr/rawCounts/genic/exonExonJunctions")
## first four columns have the exon1ID - exon2ID - chr and TSS

## Loading information of exon -exon junctions
# chr - start - end - ID
map <- read.delim(pipe("grep GB data/general/exongrps.log | cut -f2-5 -d' '"),sep=" ",header=F)
# Number of exons - length - max exon length If it's a complicated groupr and whether there are overlapping exons
mapTmp <- read.delim(pipe("grep GI data/general/exongrps.log | cut -f2-6 -d' '"),sep=" ",header=F)

map <- cbind(map,mapTmp)

## load the map of the ID with the exon IDs
mapTmp <- read.delim(pipe("grep GE data/general/exongrps.log"),sep=" ",header=F)
mapTmp[1,!is.na(mapTmp[1,])]



# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")
