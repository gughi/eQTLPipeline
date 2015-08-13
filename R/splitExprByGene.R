splitExprByGene <-
function(i,ensemblRef,ensemblGenes,PUTM,exprPUTM,SNIG,exprSNIG,pathResExpr)
{
  ## general functions ##
  sys.source("/home/adai/scripts/common_functions.R",
             attach(NULL, name="myenv"))
  
  
  ## load the ensembl file with ENSG, ENST, gene start and stop and chr
  ## ensemblRef <- read.delim(file="/home/seb/eQTL/ensemblRef.txt", as.is=T,header=T)
  
  ##load the entire list of genes
  ## ensemblGenes <- read.delim(file="/home/seb/eQTL/ensemblGenes.txt", as.is=T,header=T)
  ensemblRef.gene  <- ensemblGenes[i, ]; rm(i)
  geneID <- ensemblRef.gene$Ensembl.Gene.ID
  ## get the ensemble geneID
  gene.info <- ensemblRef[ which(ensemblRef$Ensembl.Gene.ID == geneID), ]
  
  print(geneID)
  ## "ENSG00000266195"
  
  ## filenames and folders ##
  ##pathResExpr <- "/home/seb/eQTL/expr/simpleQuantificationExonIntrons/CQNINT/13PEER/"
  dir.create(pathResExpr, showWarnings=FALSE )
  ##dir.create("/home/seb/eQTL/expr/simpleQuantificationExonIntrons/CQN/9PEER/chunks_log/", showWarnings=FALSE )
  
  fn.rda <- paste0(pathResExpr,"/byGene_snps1Mb/", geneID, ".rda")
  dir.create( dirname(fn.rda), showWarnings=FALSE )
  
  ## Read in tissues
  tissues <- c("PUTM","SNIG")
  expr <- vector(mode="list", length(tissues))
  names(expr) <- tissues
  
  ## change the sample ID to the individual ID to match the genotyped data
  #samples <- rownames(exprPUTM)
  # indID <- read.table.rows("/home/seb/phenotype/PUTMinfo.csv", keepRows=samples, sep=",")
  # IDs <- indID[samples,6]
  IDs <- PUTM[which(PUTM$A.CEL_file %in% as.character(rownames(exprPUTM))),1]
  IDs <- gsub("/","_",IDs)
  rownames(exprPUTM) <- IDs
  rm(indID,IDs)
#   samples <- rownames(exprSNIG)
#   indID <- read.table.rows("/home/seb/phenotype/SNIGinfo.csv", keepRows=samples, sep=",")
#   IDs <- indID[samples,6]
#   rownames(exprSNIG) <- IDs
  IDs <- SNIG[which(SNIG$A.CEL_file %in% as.character(rownames(exprSNIG))),1]
  IDs <- gsub("/","_",IDs)
  rownames(exprSNIG) <- IDs
  rm(IDs)
  
  
  expr[["SNIG"]] <- as.matrix(exprSNIG[,is.element(colnames(exprSNIG),ensemblRef.gene$Ensembl.Gene.ID)])
  expr[["PUTM"]] <- as.matrix(exprPUTM[,is.element(colnames(exprPUTM),ensemblRef.gene$Ensembl.Gene.ID)])
  
  ## check if the gene is expressed in one tissu or 
  if (length(expr[["SNIG"]])==0 && length(expr[["PUTM"]])==0)
  {
    cat(geneID,"\n",file=paste0(pathResExpr,"/genesNoExpressed"), append=T)
    stop("Gene not expressed")
  }
  for(tissue in tissues ){
    if (length(expr[[tissue]])==0)
    {
      cat(geneID,"\n",file=paste0(pathResExpr,"/genesNoExpressed",tissue), append=T)
      tissues <- tissues[-which(tissues==tissue)]
    }
    ## this case assign the column name to maintain the structure of the matrix, specially to keep the name of the transcript
    
    if (ncol(expr[["SNIG"]])==1)
    {
      colnames(expr[["SNIG"]]) <- intersect(colnames(exprSNIG),ensemblRef.gene$Ensembl.Gene.ID)
    }
    if (ncol(expr[["PUTM"]])==1)
    {
      colnames(expr[["PUTM"]]) <- intersect(colnames(exprPUTM),ensemblRef.gene$Ensembl.Gene.ID)
    }
  }
  
  save(expr, gene.info,
       file=fn.rda, compress="bzip2", ascii=T)
  rm(fn.rda)
  
}
