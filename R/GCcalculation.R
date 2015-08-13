GCcalculation <-
function (region,genRef) {
  ##We now define the BED file to then use it to calculate GC content
  print("creting the BED file")
  tmpBED <- tempfile("GCcont", fileext = ".BED")
  BED <- paste(gsub("chr","",region[,1]),region[,2],region[,3],rownames(region), sep="\t")
  write.table(data.frame(BED), file = tmpBED, row.names = F, 
              col.names = F, quote = F)
  rm(BED)
  
  tmpGCcon <- tempfile("GCcont")  
  cmd <- paste0("bedtools nuc -fi ",genRef," -bed ",tmpBED," > ",tmpGCcon)
  
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
