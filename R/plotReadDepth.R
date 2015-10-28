plotReadDepth <- function(gene,gen = "hg19",ensembl)
{
  ## load the library
  library(biomaRt)
  library(Gviz)
  
  ##geneType <- "lincRNA"
  
  ##gene <- "ENSG00000006555"
  
  filters <- c("chromosomal_region")
  startStop <- getBM(attributes=c("chromosome_name","start_position","end_position","external_gene_id"), filters="ensembl_gene_id", values=list(gene), mart=ensembl)
  chr <- unique(startStop$chromosome) 
  gene <- startStop$external_gene_id
  chromoReg <- paste0(chr,":",startStop$start_position,":",startStop$end_position,":-1,",
                      chr,":",startStop$start_position,":",startStop$end_position,":1")
  
  defGen <- getBM(attributes=c("chromosome_name","exon_chrom_start",
                               "exon_chrom_end","strand","gene_biotype",
                               "ensembl_gene_id",'ensembl_exon_id',
                               'ensembl_transcript_id',"external_gene_id")
                  , filters=filters, values=list(chromoReg), mart=ensembl)
  chr <- unique(defGen$chromosome) 
  rm(geneType,sta,filters)
  width <- (defGen$exon_chrom_end - defGen$exon_chrom_start) +1
  
  ## the strand field need to be converted
  ## convert 1 to + and -1 to -
  defGen <- cbind(defGen[,1:3],width,defGen[,4:9])  
  rm(width)
  
  defGen$strand <- gsub("-1","-",defGen$strand)
  defGen$strand <- gsub("1","+",defGen$strand)
  defGen$strand <- as.factor(defGen$strand)
  ## we now change the column names
  
  colnames(defGen) <- c("chromosome","start","end","width","strand",
                        "feature","gene","exon","transcript","symbol")  
  

  load(paste0("data/expr/rawCounts/intergenic/fullCoverage/fullCoverageChr",chr,".rda"))
  
  ## we load the data for the example
  
  ### Select of the bp for the gene considered
  head(fullCov)
  
  fullCovtmp <- fullCov[[1]][(startStop$start_position-10000):(startStop$end_position+10000),]
  
  rownames(fullCovtmp) <- (startStop$start_position-10000):(startStop$end_position+10000)
  
  tail(data.frame(fullCovtmp))
  
  grtrack <- GeneRegionTrack(defGen, genome = gen,
                             chromosome = chr, name = gene)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  
  ##load(file="/home/seb/")
  
  ## we load the dat
  meanCov <- apply(data.frame(fullCovtmp),1,mean)
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
  rm(meanCov)
  colnames(dat) <- c("chr","start","end","coverage")
  data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  dtrack <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",name=gene,type="histogram")
  dtrack_heat <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",name=gene,type="heatmap")
  
  ## itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  
  plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "symbol")                     
  
}
