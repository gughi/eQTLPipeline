#' ---
#' title: "eQTLs plots"
#' author: "Manuel Sebastian Guelfi"
#' ---

#' This report explains the code to visualise the eQTL on a defined region.
#' 
#' For this example we will get use the **ABLIM2** and the deletion at position **4:7983633:TG_T** and we will show it for only a brain region, in our case putamen

#+ echo=FALSE
## First we load libraries and we set the location of the data
suppressMessages(library(devtools))
setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()
suppressMessages(library(R.utils))
suppressWarnings(path <- readWindowsShortcut("data.lnk", verbose=FALSE))
setwd(dirname(path$networkPathname))
rm(path)

# First, load the gene ID from ensembl
gene <- "ENSG00000163995"
snp <- "chr4:7983633:TG_T"

# Load dosage and info about the indel
# Fisrt load all variants around the ABLIM2
load(paste0("data/snps/byGene/",gene,".rda"))
# Select the marker we are interested in
markers <- markers[snp,]

# load the samples information this is needed to select the brain region 
load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
# individuals IDs are formatted using as "_" to match the names of the file 
IDs <- gsub("/","_",PUTM$U.SD_No)
# order the dosage based on the individuals IDs
tmp <- markers[,as.character(IDs)]
names(tmp) <- as.character(PUTM$A.CEL_file)
markers <- list(info=markers[,c(1:6)],genotype=tmp)
rm(tmp,IDs)

IDs=PUTM$A.CEL_file


# load the ensembl object with version 72 to get the gene 
# definition of the gene, all exons and isoform structures
library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                   host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

## finally plot using the UKBECeQTL library
suppressWarnings(suppressMessages(
  plotLoceQTLs(gene = gene,
               ensembl = ensembl,
               IDs = IDs,
               genotype = markers)))



