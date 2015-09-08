
print( chr <- as.numeric( commandArgs(trailingOnly=T)[1] ) )

library(doParallel)
library(foreach)

cl <- makeCluster(7)
registerDoParallel(cl)
getDoParWorkers()

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
rm(covs,PUTM,RPKM.cqn)

## we load the expression data this to save the loading time

##fn <- paste0("/home/seb/imputed_v3/chr", chr, "/chr", chr, ".dose")
allMarkers <- read.delim(paste0("/home/ramasamya/genotyped/imputed_v3/chr", chr, "/chr", chr, ".dose"),sep=" ",row.names=1,check.names=FALSE)

##load the residual corrected expression
load("data/expr/normalisedCounts/intergenic/resids.PUTM.rda")

intergenicRegions <- starStopReg[which(starStopReg$chr %in% paste0("chr",chr)),]

outputFolder <- "data/snps/byRegion/PUTM/"
logFolder <- "data/snps/byRegion/logs/PUTM/"
dir.create(file.path(paste0(logFolder)), showWarnings=FALSE)


## map of the snps 
snps.map <- "/home/ramasamya/genotyped/imputed_v3/polys.map"
## path to the imputed.info file
imputed.info <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"
## path folder for regions with no SNPs    
regIDsLogNo <- "data/snps/byRegion/PUTM/tIDs_noPolys"


Sys.time()
foreach(i=1:length(intergenicRegions))%dopar%splitSNPsByRegion(i,allMarkers,intergenicRegions,outputFolder,logFolder,snps.map,imputed.info,regIDsLogNo)
Sys.time()
stopCluster(cl)


rm(outputFolder,logFolder,
   regIDsLogNo,imputed.info,
   snps.map,intergenicRegions,resids)

cl <- makeCluster(7)
registerDoParallel(cl)
getDoParWorkers()

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.SNIG")
rm(covs,SNIG,RPKM.cqn)

## we load the expression data this to save the loading time

##load the residual corrected expression
load("data/expr/normalisedCounts/intergenic/resids.SNIG.rda")

intergenicRegions <- starStopReg[which(starStopReg$chr %in% paste0("chr",chr)),]

outputFolder <- "data/snps/byRegion/SNIG/"
logFolder <- "data/snps/byRegion/logs/SNIG/"
dir.create(file.path(paste0(logFolder)), showWarnings=FALSE)


## map of the snps 
snps.map <- "/home/ramasamya/genotyped/imputed_v3/polys.map"
## path to the imputed.info file
imputed.info <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"
## path folder for regions with no SNPs    
regIDsLogNo <- "data/snps/byRegion/SNIG/tIDs_noPolys"


Sys.time()
foreach(i=1:length(intergenicRegions))%dopar%splitSNPsByRegion(i,allMarkers,intergenicRegions,outputFolder,logFolder,snps.map,imputed.info,regIDsLogNo)
Sys.time()
stopCluster(cl)





