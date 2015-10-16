
## function that return the TSS for the ensembl annotation
getTSS <- function(gene, annotation){
  
 tmp <- annotation[which(annotation$ensembl_gene_id %in% gene),] 
 if(tmp$strand <0){
    return(tmp$end_position)
 }else{
   return(tmp$start_position)
 }    
  
}