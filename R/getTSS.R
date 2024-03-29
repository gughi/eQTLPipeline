
## function that return the TSS for the ensembl annotation
getTSS <- function(gene, annotation){
  
 tmp <- annotation[which(annotation$ensembl_gene_id %in% gene),] 
 if(tmp$strand <0){
    return(tmp$end_position)
 }else{
   return(tmp$start_position)
 }    
  
}

getTES <- function(gene, annotation){
  
  tmp <- annotation[which(annotation$ensembl_gene_id %in% gene),] 
  if(tmp$strand <0){
    return(tmp$start_position)
  }else{
    return(tmp$end_position)
  }    
  
}

getDistTSS <- function(marker,gene,annotation)
{
  tmp <- annotation[which(annotation$ensembl_gene_id %in% gene),]
  if(tmp$strand <0){
    return(tmp$end_position-as.integer(marker))
  }else{
    return(as.integer(marker)-tmp$start_position)
  }
}

getDistTES <- function(marker,gene,annotation)
{
  tmp <- annotation[which(annotation$ensembl_gene_id %in% gene),]
  if(tmp$strand <0){
    return(tmp$start_position-as.integer(marker))
  }else{
    return(as.integer(marker)-tmp$end_position)
  }
}