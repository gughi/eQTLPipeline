
gene <- "ENSG00000100156"
snp <- "chr22:38457329"
##
load(paste0("data/snps/byGene/",gene,".rda"))
##
## load(paste0("/home/seb/eQTL/snps/byGene/",gene,".rda"))
markers <- markers[snp,]

# cooment 


load("data/general/sampleInfo.rda")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
IDs <- gsub("/","_",PUTM$U.SD_No)
tmp <- markers[,as.character(IDs)]
names(tmp) <- as.character(PUTM$A.CEL_file)
markers <- list(info=markers[,c(1:6)],genotype=tmp)
rm(IDs,tmp)

table(round(as.numeric(markers$genotype)))
genotype=markers
IDs=PUTM$A.CEL_file

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")


ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Feb2014.archive.ensembl.org",
                      dataset="hsapiens_gene_ensembl")


plotsplicEQTL <- function(gene,gen = "hg19",ensembl,IDs=NA, genotype)
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
  
  if(!is.na(IDs))
  {
    fullCovtmp <- fullCov[[1]][(startStop$start_position-10000):(startStop$end_position+10000),as.character(IDs)]
  }else{
    fullCovtmp <- fullCov[[1]][(startStop$start_position-10000):(startStop$end_position+10000),]
  }  
  rownames(fullCovtmp) <- (startStop$start_position-10000):(startStop$end_position+10000)
  
  tail(data.frame(fullCovtmp))
  
  
  grtrack <- GeneRegionTrack(defGen, genome = gen,
                             chromosome = chr, name = gene)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  
  ##load(file="/home/seb/")
  
  ## we load the dat

  
  tmp <- round(as.numeric(markers$genotype))
  ## we load 3expression based on teh genotype
  IDsGen <- names(markers$genotype[which(tmp %in% "0")]) 
  meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
  meanAll <- meanCov
  rm(meanCov)
  colnames(dat) <- c("chr","start","end","coverage")
  data_gRef <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  rm(dat)
  
  
  ## we load 3expression based on teh genotype
  IDsGen <- names(markers$genotype[which(tmp %in% "1")]) 
  meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
  meanAll <- cbind(meanAll,meanCov)
  rm(meanCov)
  colnames(dat) <- c("chr","start","end","coverage")
  data_gHet <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  rm(dat)
  
  ## we load 3expression based on teh genotype
  IDsGen <- names(markers$genotype[which(tmp %in% "2")]) 
  rm(tmp)
  meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
  meanAll <- cbind(meanAll,meanCov)
  rm(meanCov)
  colnames(dat) <- c("chr","start","end","coverage")
  data_gAlt <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  rm(dat)
  
  colnames(meanAll) <- c(paste0(markers$info$Al1,markers$info$Al1),
                         paste0(markers$info$Al1,markers$info$Al2),
                         paste0(markers$info$Al2,markers$info$Al2))
  


  betas <- numeric(length(rownames(fullCovtmp)))
  
  my.expr <- SlicedData$new()
  my.expr$CreateFromMatrix(as.matrix(as.data.frame(fullCovtmp[1:1000,names(genotype$genotype)])))
  my.markers <- SlicedData$new()
  my.markers$CreateFromMatrix(as.matrix(genotype$genotype))
  store <- Matrix_eQTL_main( my.markers, my.expr, output_file_name = NULL,pvOutputThreshold=1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
  head(store$all$eqtls)
    
  system.time(sapply(head(rownames(fullCovtmp),100), function(x){

    tmp <- coef(summary(glm(as.numeric(as.data.frame(fullCovtmp[x,names(genotype$genotype)])) ~ as.numeric(genotype$genotype))))
    return(tmp[2,c("Estimate","Pr(>|t|)")]) 
  }))
  
  fullCovtmp[1,names(genotype$genotype)]
  
  store$all$eqtls[which(store$all$eqtls$gene %in% "38464141"),]
  
  coef(summary(glm(as.numeric(as.data.frame(fullCovtmp[1,names(genotype$genotype)])) ~ as.numeric(genotype$genotype))))
  abline(glm(as.numeric(as.data.frame(fullCovtmp[1,names(genotype$genotype)])) ~ as.numeric(genotype$genotype)),col="red")
  
  coef(summary(glm(as.numeric(c(5,2,1)) ~ as.numeric(c(0,1,2)))))
  


      
        
  dim(as.matrix(as.data.frame(fullCovtmp[1,names(genotype$genotype)])))
  
  plot(round(as.numeric(genotype$genotype)),as.numeric(as.data.frame(fullCovtmp[1,names(genotype$genotype)])
                                                       ,xaxt="n",ylab="expression",xlab="PUTM"))
  abline(glm(PUTMtmp ~ PUTMsnp+ my.cov0[1,]+my.cov0[2,]+my.cov0[3,]),col="red")
  fit <- coef( summary(glm( ~ PUTMsnp))
  
  head(meanAll)               
               
  dtrackRef <- DataTrack(range=data_gRef,chromosome=paste0("chr",chr),genome="hg19",name=paste0(markers$info$Al1,markers$info$Al1),type="histogram",
                         ylim=c(0,max(max(data_gHet$cov),max(data_gRef$cov),max(data_gAlt$cov),na.rm=T)))
  dtrackHet <- DataTrack(range=data_gHet,chromosome=paste0("chr",chr),genome="hg19",name=paste0(markers$info$Al1,markers$info$Al2),type="histogram",
                         ylim=c(0,max(max(data_gHet$cov),max(data_gRef$cov),max(data_gAlt$cov),na.rm=T)))
  dtrackAlt <- DataTrack(range=data_gAlt,chromosome=paste0("chr",chr),genome="hg19",name=paste0(markers$info$Al2,markers$info$Al2),type="histogram",
                         ylim=c(0,max(max(data_gHet$cov),max(data_gRef$cov),max(data_gAlt$cov),na.rm=T)))
  
  rm(data_gHet,data_gAlt,data_gRef)
  ## itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  class(dtrackRef)
  tracks <- list(gtrack)
  if(nrow(values(dtrackRef))>0)
  {
    tracks <- list(c(unlist(tracks),dtrackRef))
  }
  if(nrow(values(dtrackHet))>0)
  {
    tracks <- list(c(unlist(tracks),dtrackHet))
  }
  if(nrow(values(dtrackAlt))>0)
  {
    tracks <- list(c(unlist(tracks),dtrackAlt))
  }
  tracks <- list(c(unlist(tracks),grtrack))
  plotTracks(unlist(tracks),transcriptAnnotation = "symbol")
  
  colnames(meanAll) <- c(paste0(markers$info$Al1,markers$info$Al1),paste0(markers$info$Al1,markers$info$Al2),paste0(markers$info$Al2,markers$info$Al2))
  pvalstrack <- DataTrack(data=t(meanAll[,2:3]),start=as.numeric(rownames(meanAll)),end=as.numeric(rownames(meanAll)), chromosome = chr, genome = gen,
                          name = "Both",type=c("p"),groups=c("GA","AA"),col=c("black","red"))
  
  plotTracks(list(gtrack, pvalstrack,grtrack),transcriptAnnotation = "symbol",type="a")  
  
  ## The code below needs to improve
  #   colnames(meanAll) <- c(paste0(markers$info$Al1,markers$info$Al1),paste0(markers$info$Al1,markers$info$Al2),paste0(markers$info$Al2,markers$info$Al2))
  #   allInSameTrack <- DataTrack(data=t(meanAll),start=as.numeric(rownames(meanAll)),end=as.numeric(rownames(meanAll)), chromosome = chr, genome = gen,
  #                               name = "All",type=c("p"),groups=c(paste0(markers$info$Al1,markers$info$Al1),
  #                                                                 paste0(markers$info$Al1,markers$info$Al2),
  #                                                                 paste0(markers$info$Al2,markers$info$Al2))
  #                               ,col=c("black","red","blue"))
  #   
  #   plotTracks(list(gtrack, allInSameTrack,grtrack),transcriptAnnotation = "symbol",type="a")  
  
  
}
