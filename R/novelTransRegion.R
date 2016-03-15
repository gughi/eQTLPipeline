
## This functions identify the novel transcribed regions in intragenic
novelTransRegion <- function(geneInfo,ensembl,threshold=10,tissue="PUTM")
{
  load(paste0("data/expr/rawCounts/intergenic/transcribedRegions/",tissue,"/chr",geneInfo[3],".rda"))
  rm(annotatedRegions,filteredCov,intergenicRegions,coverageIntergenic)
  library(GenomicRanges)
  library(biomaRt)
  ## we subset the regions getting the regions overlapping to the gene 
  #     subsetByOverlaps(expressedRegions$chr1$regions, GRanges(paste0("chr",geneInfo[3]),
  #                                                             IRanges(as.numeric(geneInfo[4]), as.numeric(geneInfo[5]))))
  
  #head(expressedRegions)
  
  
  ## select the neuro genes that non overlap any other gene
  
  exonDef <- getBM(attributes=c("ensembl_gene_id","chromosome_name","ensembl_exon_id","exon_chrom_start","exon_chrom_end"),
                   filters="ensembl_gene_id",
                   values=geneInfo[1], mart=ensembl)
  
  
  exonDef <- GRanges(paste0("chr",exonDef[,2]), IRanges(exonDef[,4], exonDef[,5]))
  ## select the transcribed regions identified in the data
  tmp <- subsetByOverlaps(expressedRegions[[1]]$regions,
                          GRanges(paste0("chr",geneInfo[3]),
                                  IRanges(as.numeric(geneInfo[4]), as.numeric(geneInfo[5]))))
  
  
  tmp2 <- countOverlaps(tmp[which(tmp$value>threshold),], exonDef)==0
  tmp2 <- tmp2[which(tmp2 %in% TRUE)]
    
  ##rm(expressedRegions)
  return(list(gene=geneInfo[1],
          totalCovered=length(tmp),
          coveredNoAnn=table(countOverlaps(tmp[which(tmp$value>threshold),], exonDef)==0)["TRUE"],
          expressedRegions[[1]]$regions[names(tmp2)],
          totalExons=length(reduce(exonDef))))
  
}


## function that returns, if any, the unannotated regions for the gene of interest 
getUnAnnotatedRegions <- function(gene,ensembl,tissue=NA,threshold=10)
{
  
  library(biomaRt)
  exonDef <- getBM(attributes=c("ensembl_gene_id","chromosome_name","ensembl_exon_id","exon_chrom_start","exon_chrom_end"),
                   filters="ensembl_gene_id",
                   values=gene, mart=ensembl)
  
  ## we load the results from derfinder 
  load(paste0("data/expr/rawCounts/intergenic/transcribedRegions/",tissue,"/chr",unique(exonDef$chromosome_name),".rda"))
  rm(annotatedRegions,filteredCov,intergenicRegions,coverageIntergenic)
  library(GenomicRanges)
  
  
  ## select only the coverage for the genes
  geneDef <- GRanges(paste0("chr",unique(exonDef[,2])), IRanges(min(exonDef[,4]), max(exonDef[,5])))  
  
  ## select the transcribed regions identified in the data
  tmp <- subsetByOverlaps(expressedRegions[[1]]$regions,geneDef)
  rm(expressedRegions)
  tmp <- tmp[which(tmp$value>threshold),]
  
  
  exonDef <- GRanges(paste0("chr",exonDef[,2]), IRanges(exonDef[,4], exonDef[,5]))
  
  tmp2 <- countOverlaps(tmp[which(tmp$value>10),], exonDef)==0
  tmp2 <- tmp2[which(tmp2 %in% TRUE)]
  
  return(tmp[names(tmp2),])
  
}

