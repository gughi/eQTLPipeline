
plotReadDepth <- function(gene,gen = "hg19",ensembl,IDs=NA,chr=NA,from=NULL,to=NULL)
{
  ## load the library
  library(biomaRt)
  library(Gviz)
  library(GenomicRanges)
  
  ##geneType <- "lincRNA"
  
  ##gene <- "ENSG00000006555"
  
  filters <- c("chromosomal_region")
  if (!is.na(gene))
  {
    
    startStop <- getBM(attributes=c("chromosome_name","start_position","end_position","external_gene_id"), filters="ensembl_gene_id", values=list(gene), mart=ensembl)
    chr <- unique(startStop$chromosome) 
    gene <- startStop$external_gene_id
    
    
  }else{
    startStop <- as.data.frame(t(as.matrix(c(chr,from,to,""))))
    colnames(startStop) <- c("chromosome_name","start_position","end_position","external_gene_id")
    startStop$start_position <- as.numeric(as.character(startStop$start_position))
    startStop$end_position <- as.numeric(as.character(startStop$end_position))
    print(startStop)
  }
  
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
  


  grtrack <- GeneRegionTrack(defGen, genome = gen,
                             chromosome = paste0("chr",chr), name = as.character(startStop$external_gene_id))
  gtrack <- GenomeAxisTrack()
  ##itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  
  ##load(file="/home/seb/")
  
  ## we load the dat
  meanCov <- apply(data.frame(fullCovtmp),1,mean)
  meanCov[which(meanCov < 10)]=0 
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))  
  rm(meanCov)
  print(head(dat))
  colnames(dat) <- c("chr","start","end","coverage")
  
  data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  #print(head(data_g))
  dtrack <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",name="coverage",type="histogram")
  #dtrack <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",name="SSS",type="histogram")
  #print(dtrack)
  dtrack_heat <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",name=as.character(startStop$external_gene_id),type="heatmap")
  
  ## itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  ##print(head(dtrack))
  
  #plotTracks(list(gtrack,dtrack),transcriptAnnotation = "symbol")
  plotTracks(list(gtrack, dtrack,grtrack,dtrack_heat),transcriptAnnotation = "symbol")      
  
  
  
#   dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
#   rm(meanCov)
#   colnames(dat) <- c("chr","start","end","coverage")
#   data_g <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
#   dtrack <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",type="histogram",name="coverage")
#   dtrack_heat <- DataTrack(range=data_g,chromosome=paste0("chr",chr),genome="hg19",type="heatmap")
  
  
  
}

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
  #itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  
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



plotLoceQTLs <- function(gene,gen = "hg19",ensembl,IDs=NA, genotype, highLight=NULL,from=NULL,to=NULL)
{
  ## load the library
  library(biomaRt)
  library(Gviz)
  
  ##geneType <- "lincRNA"
  
  ##gene <- "ENSG00000006555"
  
  filters <- c("chromosomal_region")
  if (!is.na(gene))
  {
    
    startStop <- getBM(attributes=c("chromosome_name","start_position","end_position","external_gene_id"), filters="ensembl_gene_id", values=list(gene), mart=ensembl)
    chr <- unique(startStop$chromosome) 
    gene <- startStop$external_gene_id
    
  
  }else{
    chr <- gsub("chr","",unlist(strsplit(rownames(genotype$info),":"))[1])
    startStop <- as.data.frame(t(as.matrix(c(chr,from,to,""))))
    colnames(startStop) <- c("chromosome_name","start_position","end_position","external_gene_id")
    startStop$start_position <- as.numeric(as.character(startStop$start_position))
    startStop$end_position <- as.numeric(as.character(startStop$end_position))
    print(startStop)
  }
  
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
                             chromosome = chr, name = as.character(startStop$external_gene_id))
  gtrack <- GenomeAxisTrack()
  ##itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  
  ##load(file="/home/seb/")
  
  
  ## we load the dat
  
  
  tmp <- round(as.numeric(genotype$genotype))
  ## we load 3expression based on teh genotype
  IDsGen <- names(genotype$genotype[which(tmp %in% "0")]) 
  meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
  meanAll <- meanCov
  rm(meanCov)
  colnames(dat) <- c("chr","start","end","coverage")
  data_gRef <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  rm(dat)
  
  
  ## we load 3expression based on teh genotype
  IDsGen <- names(genotype$genotype[which(tmp %in% "1")]) 
  meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
  meanAll <- cbind(meanAll,meanCov)
  rm(meanCov)
  colnames(dat) <- c("chr","start","end","coverage")
  data_gHet <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  rm(dat)
  
  ## we load 3expression based on teh genotype
  IDsGen <- names(genotype$genotype[which(tmp %in% "2")]) 
  rm(tmp)
  meanCov <- apply(data.frame(fullCovtmp[,IDsGen]),1,mean)
  dat <- data.frame(paste0("chr",chr),as.numeric(rownames(fullCovtmp)),as.numeric(rownames(fullCovtmp)),as.vector(meanCov))                     
  ## filter
  #meanCov[which(meanCov < 10)]=0 
  
  meanAll <- cbind(meanAll,meanCov)
  rm(meanCov)
  colnames(dat) <- c("chr","start","end","coverage")
  data_gAlt <- with(dat, GRanges(chr, IRanges(start, end), cov=coverage))
  rm(dat)
  
  
  colnames(meanAll) <- c(paste0(genotype$info$Al2,genotype$info$Al2),
                         paste0(genotype$info$Al1,genotype$info$Al2),
                         paste0(genotype$info$Al1,genotype$info$Al1))
  print(head(meanAll))
  
  library(MatrixEQTL)
  my.expr <- SlicedData$new()
  my.expr$CreateFromMatrix(as.matrix(as.data.frame(fullCovtmp[,names(genotype$genotype)])))
  my.markers <- SlicedData$new()
  my.markers$CreateFromMatrix(as.matrix(genotype$genotype))
  store <- Matrix_eQTL_main( my.markers, my.expr, output_file_name = NULL,pvOutputThreshold=1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
  head(store$all$eqtls)
  
  ## betas
  datBetas <- GRanges(chr, IRanges(as.numeric(as.character(store$all$eqtls$gene))
                                   , as.numeric(as.character(store$all$eqtls$gene))), cov=store$all$eqtls$beta)
  
  dtrackBetas <- DataTrack(range=datBetas,chromosome=paste0("chr",chr),genome="hg19",name="betas",type="histogram")
  
  datPval <- GRanges(chr, IRanges(as.numeric(as.character(store$all$eqtls$gene))
                                  , as.numeric(as.character(store$all$eqtls$gene))), cov=-log10(store$all$eqtls$pvalue))
  
  dtrackPval <- DataTrack(range=datPval,chromosome=paste0("chr",chr),genome="hg19",name="-log10(pvalues)",type="gradient")                          
  
  ## converted to get the right oreder when displaying the layers 
  meanAll <- t(meanAll)[c(paste0(genotype$info$Al1,genotype$info$Al1),
                          paste0(genotype$info$Al1,genotype$info$Al2),
                          paste0(genotype$info$Al2,genotype$info$Al2)),]
  
  ## The code below needs to improve
  ##colnames(meanAll) <- c(paste0(genotype$info$Al2,genotype$info$Al2),paste0(genotype$info$Al1,genotype$info$Al2),paste0(genotype$info$Al1,genotype$info$Al1))
  allInSameTrack <- DataTrack(data=meanAll,start=as.numeric(colnames(meanAll)),end=as.numeric(colnames(meanAll)), chromosome = chr, genome = gen,
                              name = "Stratified raw counts",type=c("l"),groups=rownames(meanAll),legend = TRUE
                              ,col=c("black","red","blue"))
  #
  if(!is.null(highLight)){
    #        ht <- HighlightTrack(trackList = list(allInSameTrack,grtrack), 
    #                             start = highLight[,2], 
    #                             width = (highLight[,3]-highLight[,2]),chromosome = as.numeric(highLight[,1]),
    #                             col=scales::alpha(c(1:nrow(highLight)),.5),inBackground=TRUE,
    #                             fill=scales::alpha(c(1:nrow(highLight)),.5))
    #     
    
    ## This code highlight the exon exon junctions 
    
    
    highLight <- rbind(highLight[,c("chrExon1","startExon1","endExon1","exExID")],
                       setNames( highLight[,c("chrExon1","startExon2","endExon2","exExID")],
                                 c("chrExon1","startExon1","endExon1","exExID")))
    
    highLight$width <- highLight[,3]-highLight[,2] 
    highLight$symbol <- gene 
    highLight$strand <- "*"
    highLight <- highLight[,c("chrExon1","startExon1","endExon1","strand","width","exExID","symbol")]    
    colnames(highLight) <- c("chromosome","start","end","strand","width","transcript", "symbol")
    
    
    exExJuntrack <- GeneRegionTrack(highLight, genome = gen,
                                    chromosome = chr, name = "exon-exon junction",col="red",fill="red")
    
    
    plotTracks(list(gtrack,dtrackPval,dtrackBetas,allInSameTrack,grtrack,exExJuntrack),transcriptAnnotation = "transcript")
    
  }else{
    
    plotTracks(list(gtrack,dtrackPval,dtrackBetas,allInSameTrack,grtrack),transcriptAnnotation = "symbol",from = 23322001 , to=23322800)
    #plotTracks(list(gtrack,dtrackPval,dtrackBetas,allInSameTrack,grtrack),transcriptAnnotation = "symbol",from = from, to=to)
  }
  
  
}





