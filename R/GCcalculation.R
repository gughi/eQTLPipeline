#' Calculate the GC content
#' 
#' @param data frame with the regions with the following fixed structure: first column is the chr, second column is the start position, third column is the stop position, forth column is the Identifier
#' @param the location of the .fa file of the genome reference to calculate the GC content
#' @note this function need bedtools installed in your system
#' @return data frame with one column as GC content, the rownames are the identifiers
GCcalculation <-
function (region,genRef,pathBedtools) {
  ##We now define the BED file to then use it to calculate GC content
  print("creting the BED file")
  tmpBED <- tempfile("GCcont", fileext = ".BED")
  BED <- paste(gsub("chr","",region[,1]),region[,2],region[,3],rownames(region), sep="\t")
  write.table(data.frame(BED), file = tmpBED, row.names = F, 
              col.names = F, quote = F)
  rm(BED)
  
  tmpGCcon <- tempfile("GCcont")  
  cmd <- paste0(pathBedtools," nuc -fi ",genRef," -bed ",tmpBED," > ",tmpGCcon)
  
  print("executing the bedtools")
  system(cmd)
  rm(cmd)
  
  print("collecting the GC content")
  GCcontent <- read.delim(pipe(paste("cut -f4,6", tmpGCcon)))
  colnames(GCcontent) <- c("ID","GCcontent")
  rownames(GCcontent) <- GCcontent$ID
  GCcontent$ID <- NULL 
  GCcontent
}
