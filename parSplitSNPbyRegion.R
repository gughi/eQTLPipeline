chr<- 21


library(doParallel)
library(foreach)

detectCores()
## [1] 24

cl <- makeCluster(16)
registerDoParallel(cl)
getDoParWorkers()

load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")
rm(covs,PUTM,RPKM.cqn)

## we load the expression data this to save the loading time

##fn <- paste0("/home/seb/imputed_v3/chr", chr, "/chr", chr, ".dose")
allMarkers <- read.delim(paste0("/home/ramasamya/genotyped/imputed_v3/chr", chr, "/chr", chr, ".dose"),sep=" ",row.names=1,check.names=FALSE)

##load the residual corrected expression
load("data/expr/normalisedCounts/intergenic/resids.PUTM.rda")

outputFolder <- "data/snps/byRegion/PUTM/"
logFolder <- "data/snps/byRegion/logs/PUTM/"
dir.create(paste0(outputFolder,"chr",chr), showWarnings=FALSE )
dir.create(file.path(paste0(logFolder)), showWarnings=FALSE)


## map of the snps 
snps.map <- "/home/ramasamya/genotyped/imputed_v3/polys.map"
## path to the imputed.info file
imputed.info <- "/home/ramasamya/genotyped/imputed_v3/imputed.info"
## path folder for regions with no SNPs    
regIDsLogNo <- "data/snps/byRegion/PUTM/tIDs_noPolys/"


Sys.time()
foreach(i=1:length(intergenicRegions))%dopar%splitSNPsByRegion(i,chr,allMarkers,intergenicRegions,outputFolder,snps.map,imputed.info,regIDsLogNo)
Sys.time()
stopCluster(cl)







