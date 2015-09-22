

## load the residual corrected expression
load("data/expr/normalisedCounts/genic/exons/resids.PUTM.rda")
## load the information about tissue
load("data/expr/normalisedCounts/genic/exons/RPKM.cqn.PUTM")
## remove the object we don't need
rm(covs,RPKM.cqn)


load("data/general/sampleInfo.rda")
IDs <- sampleInfo[which(sampleInfo$A.CEL_file %in% as.character(rownames(resids))),"U.SD_No"]
IDs <- gsub("/","_",IDs)
rownames(resids) <- IDs
rm(indID,IDs)

snpLocation <- "/home/guelfi/eQTL/snps/byGene/"

outputFolder <- "/home/guelfi/eQTLPipeline/data/results/genic/exons/resMatrixEQTL/PUTM/"
fullResults <- "/home/guelfi/eQTLPipeline/data/results/genic/exons/fullResults/PUTM/"

genotypeFile <- "/home/guelfi/plinkOutput/eigenvec"


my.covTMP <- read.table.rows(paste0("/home/guelfi/plinkOutput/eigenvec"), keepRows=rownames(resids), sep=" ",header=F)

Sys.time()
foreach(i=1:ncol(resids))%dopar%runCisEQTLExons(i=i,resids=resids,
                                                     snpLocation=snpLocation,
                                                     outputFolder=outputFolder,
                                                     my.covTMP=my.covTMP,
                                                     fullResults=fullResults)

Sys.time()
stopCluster(cl)
rm(cl)

