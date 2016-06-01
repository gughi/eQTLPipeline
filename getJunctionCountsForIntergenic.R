## scripts that checks whether the intergenic regions identified has any exon-exon junction connected to the nearest gene.
## it returns a matrix with boolean value that tells whether it has or not the exon-exon junction and the counts


library(devtools)
##setwd("C:/Users/mguelfi/projectsR/eQTLPipeline/")
load_all()


library(R.utils)
path <- readWindowsShortcut("data.lnk", verbose=FALSE)
setwd(dirname(path$networkPathname))
rm(path)

## load teh intergenic regions
load(file="data/general/defIntergenic.rda")

## we remove Intergenic regions that don't have any new portion but that are totally inside the gene
tmp <- regInformation[which(regInformation$uniquePortion > 0),]
tmp$overlap <- as.numeric(as.character(tmp$overlap))
tmp <- tmp[which(tmp$overlap < 2),]
rm(regInformation)

genesMap <- read.delim("data/general/ensemblGenes.txt")
## define the genomic regions
GR <- GRanges(seqnames = Rle(genesMap$Chromosome.Name),
              ranges = IRanges(start=genesMap$Gene.Start..bp.,
                               end = genesMap$Gene.End..bp.,
                               names = genesMap$Ensembl.Gene.ID))


bedPath <- list.files("data/junctions/PUTM/")

print(paste("Start Time",Sys.time()))
for (i in 1:length(bedPath))
{
  intNovelExon <- identifyNewExon(interTable = tmp,juncPath = paste0("data/junctions/PUTM/",bedPath[i]),GR = GR)
  save(intNovelExon,file=paste0(as.character("data/junctions/PUTM/"),unlist(strsplit(as.character(bedPath[i]),"junctions"))[1],"intergenicNovelExon.rda"))
  print(paste0(i/length(bedPath)*100,"% files Completed ",Sys.time()))
}











    

