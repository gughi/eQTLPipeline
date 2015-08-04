## This script contains the general gunction needed for the eQTL analysis

## calculate the GCcontent
expr <- expr[which(expr$width >= 100),]

## now we restricted our regions to 45114

## We need to calculate the GC content for these regions

### for trackable purposes we will give ID to regions the ID would be DER+number 
IDs <- paste0("DER",c(1:nrow(expr)))
rownames(expr) <- IDs
rm(IDs)

##We now define the BED file to then use it to calculate GC content
BED <- paste(gsub("chr","",expr$chr),expr$start,expr$end,rownames(expr), sep="\t")

write.table(data.frame(BED), file = paste0("/home/seb/DERFINDER/eQTL/PUTMSNIGGCcontent.BED"), row.names = F, 
            col.names = F, quote = F)
rm(BED)

cmd <- paste0("bedtools nuc -fi /home/seb/reference/genome37.72.fa -bed /home/seb/DERFINDER/eQTL/PUTMSNIGGCcontent.BED > /home/seb/DERFINDER/eQTL/PUTMSNIGGCcontent")
system(cmd)
rm(cmd)
## remove tmp files
system("rm /home/seb/DERFINDER/eQTL/PUTMSNIGGCcontent.BED")

GCcontent <- read.delim(pipe("cut -f4,6 /home/seb/DERFINDER/eQTL/PUTMSNIGGCcontent"))
colnames(GCcontent) <- c("ID","GC")
rownames(GCcontent) <- GCcontent$ID
GCcontent$ID <- NULL 
