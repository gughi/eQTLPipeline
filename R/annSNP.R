


annSNP <- function(eQTLres)
{
    library(biomaRt)
    ensembl <- useMart(biomart="ENSEMBL_MART_SNP", host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_snp")

    snps <- NULL
    snps$chr <- gsub("chr","",unlist(lapply(strsplit(as.character(eQTLres$snps),":")
                                            ,function(x){c(x[1])})))
    snps$start <- unlist(lapply(strsplit(as.character(eQTLres$snps),":")
                                ,function(x){c(x[2])}))
    
    snps$end <- unlist(lapply(strsplit(as.character(eQTLres$snps),":")
                              ,function(x){c((as.numeric(x[2])+2))}))

    snps <- as.data.frame(snps)

    rs <- apply(snps,1,function(x){
      rsID <- getBM(attributes = c("refsnp_id", "allele","chr_name","chrom_start"),
                    filters = c("chr_name", "chrom_start", "chrom_end"),
                    values=list(x[1],x[2],x[3]),
                    mart = ensembl)
      marPos <- paste0("chr",x[1],":",x[2])
      rsID <- cbind(rsID,marPos)
      return(rsID)
    }
    )

    library (plyr)
    df <- ldply (rs, data.frame)

    rs <- sapply(eQTLres$snps,function(x){
      snpPos <- paste(unlist(strsplit(as.character(x),":"))[1:2],collapse=":")
      tmp <- df[which(df$marPos %in%snpPos),]
      c(paste(tmp$refsnp_id,collapse = ";"),paste(tmp$allele,collapse = ";"))
    })

    eQTLres <- cbind(eQTLres,t(rs))
    rm(df,rs)    
    return(eQTLres)
}