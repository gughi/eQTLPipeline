## summary of eQTLs

###################
### gene exonic ###
###################

load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")

load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")

numeQTLs10 <- c(nrow(eQTLPUTM),nrow(eQTLSNIG))
names(numeQTLs10) <- c("PUTM","SNIG")
barplot(numeQTLs10,col='skyblue',main= "eQTLs gene exonic",border=F,
        sub=paste("PUTM:",numeQTLs10["PUTM"],"SNIG:",numeQTLs10["SNIG"]))

numeQTLs5 <- c(length(which(eQTLPUTM$myFDR<0.05)),length(which(eQTLSNIG$myFDR<0.05)))
names(numeQTLs5) <- c("PUTM","SNIG")
barplot(numeQTLs5,add=T,col=scales::alpha('red',.5),border=F)


numeQTLs1 <- c(length(which(eQTLPUTM$myFDR<0.01)),length(which(eQTLSNIG$myFDR<0.01)))
names(numeQTLs1) <- c("PUTM","SNIG")
barplot(numeQTLs1,add=T,col=scales::alpha('green',.5),border=F)
legend("topright",c("10%","5%","1%"),col=c('skyblue','red','green'),title="FDR",pch=15)

## table of the independecy of the signals
sign <- cbind(table(eQTLPUTM$degree),table(eQTLSNIG$degree))

par(mar=c(9,3,1,1))
barplot(sort(table(eQTLPUTM$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLPUTM$gene_biotype))),las=2)

barplot(sort(table(eQTLSNIG$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLSNIG$gene_biotype))),las=2,main="test")

rm(numeQTLs1,numeQTLs10,numeQTLs5)

############################
### gene exonic-Intronic ###
############################


load("data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")

load("data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")


numeQTLs10 <- c(nrow(eQTLPUTM),nrow(eQTLSNIG))
names(numeQTLs10) <- c("PUTM","SNIG")
barplot(numeQTLs10,col='skyblue',main= "eQTLs gene exonic-intronic",,border=F,
        sub=paste("PUTM:",numeQTLs10["PUTM"],"SNIG:",numeQTLs10["SNIG"]))

numeQTLs5 <- c(length(which(eQTLPUTM$myFDR<0.05)),length(which(eQTLSNIG$myFDR<0.05)))
names(numeQTLs5) <- c("PUTM","SNIG")
barplot(numeQTLs5,add=T,col=scales::alpha('red',.5),border=F)

numeQTLs1 <- c(length(which(eQTLPUTM$myFDR<0.01)),length(which(eQTLSNIG$myFDR<0.01)))
names(numeQTLs1) <- c("PUTM","SNIG")
barplot(numeQTLs1,add=T,col=scales::alpha('green',.5),border=F)
legend("topright",c("10%","5%","1%"),col=c('skyblue','red','green'),title="FDR",pch=15)

sign <- cbind(table(eQTLPUTM$degree),table(eQTLSNIG$degree))

sign
rm(numeQTLs1,numeQTLs10,numeQTLs5)


##################
### intergenic ###
##################


load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")

load("data/results/finaleQTLs/intergenic.Ann.SNIG.rda")


numeQTLs10 <- c(nrow(eQTLPUTM),nrow(eQTLSNIG))
names(numeQTLs10) <- c("PUTM","SNIG")
barplot(numeQTLs10,col='skyblue',main= "eQTLs intergenic",border=F,
        sub=paste("PUTM:",numeQTLs10["PUTM"],"SNIG:",numeQTLs10["SNIG"]))

numeQTLs5 <- c(length(which(eQTLPUTM$myFDR<0.05)),length(which(eQTLSNIG$myFDR<0.05)))
names(numeQTLs5) <- c("PUTM","SNIG")
barplot(numeQTLs5,add=T,col=scales::alpha('red',.5),border=F)

numeQTLs1 <- c(length(which(eQTLPUTM$myFDR<0.01)),length(which(eQTLSNIG$myFDR<0.01)))
names(numeQTLs1) <- c("PUTM","SNIG")
barplot(numeQTLs1,add=T,col=scales::alpha('green',.5),border=F)
legend("topright",c("10%","5%","1%"),col=c('skyblue','red','green'),title="FDR",pch=15)

sign <- cbind(table(eQTLPUTM$degree),table(eQTLSNIG$degree))

rm(numeQTLs1,numeQTLs10,numeQTLs5)



##################################
## Q-Q plots of random pvalues ###
##################################


library(devtools)
load_all()

## load teh random pvalues gene exonic 
load("data/results/finaleQTLs/randomPvals/pValGeneExo.PUTM.rda")
jpeg("plots/Q-QPUTMGeneeoxnic.jpeg")
qq.plot(pvalues,main="Q-Q plot of pvalues PUTM gene Exonic")
dev.off()
load("data/results/finaleQTLs/randomPvals/pValGeneExo.SNIG.rda")
jpeg("plots/Q-QSNIGGeneeoxnic.jpeg")
qq.plot(pvalues,main="Q-Q plot of pvalues SNIG gene Exonic")
dev.off()

## load teh random pvalues gene intornic+exonic 
load("data/results/finaleQTLs/randomPvals/pValGeneExoIntron.PUTM.rda")
jpeg("plots/Q-QPUTMGeneExonicIntronic.jpeg")
qq.plot(pvalues,main="Q-Q plot of pvalues PUTM gene Exonic")
dev.off()
load("data/results/finaleQTLs/randomPvals/pValGeneExo.SNIG.rda")
jpeg("plots/Q-QSNIGGeneExonicIntronic.jpeg")
qq.plot(pvalues,main="Q-Q plot of pvalues SNIG gene Exonic")
dev.off()


## load teh random pvalues intergenic
load("data/results/finaleQTLs/randomPvals/pValIntergenic.PUTM.rda")
jpeg("plots/Q-QPUTMIntergenic.jpeg")
qq.plot(pvalues,main="Q-Q plot of pvalues PUTM Intergenic")
dev.off()
load("data/results/finaleQTLs/randomPvals/pValIntergenic.SNIG.rda")
jpeg("plots/Q-QSNIGIntergenic.jpeg")
qq.plot(pvalues,main="Q-Q plot of pvalues SNIG gene Intergenic")
dev.off()




###############################
## look in the gwas overlap ###
###############################

library (plyr)


gwasCat <- read.delim("data/general/GWAS_Disease_Classification_070415.tsv")
my.eQTLs   <- read.delim(file="data/results/finaleQTLs/eQTLintergenic.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(my.eQTLs$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
head(df)
rownames(df) <- NULL
my.eQTLs <- df[which(df[,7] < 0.01),]
dim(my.eQTLs)
hits <- gsub("chr","",my.eQTLs[,1])

cat(hits, sep="\n", file=(tmpf <- tempfile()))
map <- read.delim( pipe(paste0("fgrep -w --file=", tmpf, " /home/ramasamya/genotyped/imputed_v3/chrpos2rsid_UKBEC")), sep=" ", header=FALSE, as.is=TRUE )
hits <- map$V2[!is.na(map$V2)]

hitsGwas <- unlist(lapply(strsplit(as.character(gwasCat$Strongest.SNP.Risk.Allele),"-"),function(x){x[1]}))

out <- gwasCat[which(hitsGwas %in% hits),]

print(nrow(out))
print(paste("number of SNPs found in the GWAS catalogue related with brain:",nrow(out[which(out$General.Phenotype=="brain"),])))
print(paste("number of SNPs found in the GWAS catalogue related with Adult neurological disease:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])))
print(paste("rate Adult Neuro Br:",nrow(out[which(out$General.Phenotype=="brain"& out$Adult.Neurological.disorder=="yes"),])/nrow(out)))
print(paste("rate:",nrow(out[which(out$General.Phenotype=="brain"),])/nrow(out)))


#########################################
#### Looked at the disease ontology  ####
#########################################


library(DOSE)

rm(list=ls())

load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.01),]

library("biomaRt")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

entrezID <- getBM(attributes=c("entrezgene"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=eQTLPUTM$gene, mart=ensembl)

DOres <- enrichDO(entrezID[,1], pvalueCutoff=0.01)
head(summary(DOres))



entrezID <- getBM(attributes=c("entrezgene"),
                  verbose = T,
                  filters="ensembl_gene_id",
                  values=eQTLSNIG$gene, mart=ensembl)


DOres <- enrichDO(entrezID[,1], pvalueCutoff=0.01)
head(summary(DOres))

## intronic qunatification



#############################################
##### share between diff qunatifications ####
#############################################

load(file="data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
load(file="data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")

eQTLSNIGEx <- eQTLSNIG
eQTLPUTMEx <- eQTLPUTM

load(file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
load(file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")

library(gplots)
## FDR10%
venn(list(PUTMEI=unique(eQTLPUTM$gene),SNIGEI=unique(eQTLSNIG$gene),SNIGE=unique(eQTLSNIGEx$gene),PUTME=unique(eQTLPUTMEx$gene)))
venn(list(PUTMEI=unique(eQTLPUTM$gene),PUTME=unique(eQTLPUTMEx$gene)))
venn(list(SNIGEI=unique(eQTLSNIG$gene),SNIGE=unique(eQTLSNIGEx$gene)))


## FDR5%    
eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.05),]
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.05),]

eQTLSNIGEx <- eQTLSNIGEx[which(eQTLSNIGEx$myFDR<0.05),]
eQTLPUTMEx <- eQTLPUTMEx[which(eQTLPUTMEx$myFDR<0.05),]


venn(list(PUTMEI=unique(eQTLPUTM$gene),SNIGEI=unique(eQTLSNIG$gene),SNIGE=unique(eQTLSNIGEx$gene),PUTME=unique(eQTLPUTMEx$gene)))
venn(list(PUTMEI=unique(eQTLPUTM$gene),PUTME=unique(eQTLPUTMEx$gene)))
venn(list(SNIGEI=unique(eQTLSNIG$gene),SNIGE=unique(eQTLSNIGEx$gene)))



## FDR1%    
eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.01),]
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]

eQTLSNIGEx <- eQTLSNIGEx[which(eQTLSNIGEx$myFDR<0.01),]
eQTLPUTMEx <- eQTLPUTMEx[which(eQTLPUTMEx$myFDR<0.01),]


venn(list(PUTMEI=unique(eQTLPUTM$gene),SNIGEI=unique(eQTLSNIG$gene),SNIGE=unique(eQTLSNIGEx$gene),PUTME=unique(eQTLPUTMEx$gene)))

venn(list(PUTMEI=unique(eQTLPUTM$gene),PUTME=unique(eQTLPUTMEx$gene)))
venn(list(SNIGEI=unique(eQTLSNIG$gene),SNIGE=unique(eQTLSNIGEx$gene)))


##########################################
### Relationship between the eQTLs PUTM ##
##########################################


load("data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
eQTLPUTMEI <- eQTLPUTM
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTMI <- eQTLPUTM
load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")
eQTLPUTMInter <- eQTLPUTM
eQTLPUTMExons <- read.delim("data/results/finaleQTLs/exons.PUTM.txt")
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTMEx <- eQTLPUTM
rm(eQTLPUTM)

library("gplots")
venn(list(PUTMEI=unique(eQTLPUTMEI$gene),
          PUTMEx=unique(eQTLPUTMEx$gene),
          PUTMIntr=unique(eQTLPUTMI$gene)))

## overlap gene Exons and intronic

venn(list(PUTMEx=unique(eQTLPUTMEx$gene),
          PUTMIntr=unique(eQTLPUTMI$gene)))

## overlap gene Exons+Intronic and intronic 

venn(list(PUTMEI=unique(eQTLPUTMEI$gene),
          PUTMIntr=unique(eQTLPUTMI$gene)))


plot(density(eQTLPUTMI[which(eQTLPUTMI$myFDR<0.01),"beta"]),col='red',main = "Density beta values PUTM")
lines(density(eQTLPUTMEx[which(eQTLPUTMEx$myFDR<0.01),"beta"]),col='skyblue')
legend("topright",c("Whole-gene Exonic","Whole-gene Intronic"),fill=c("skyblue","red"))


plot(density(eQTLPUTMEI[which(eQTLPUTMEI$myFDR<0.01),"beta"]))



##################################
#### Comparison with Lappainen ###
##################################



load(file="data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
load(file="data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")


## FDR5%    
eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.05),]
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.05),]


cmd <- "cut -f1,3 /home/guelfi/eQTL/465lymphoblastoid/EUR373.gene.cis.FDR5.all.rs137.txt | cut -d '.' -f1 | cut -f1,2"
lym.eQTLs.gene <- read.delim(pipe(cmd),header=T)
lym.eQTLs <- read.delim(file="/home/guelfi/eQTL/465lymphoblastoid/EUR373.gene.cis.FDR5.all.rs137.txt",header=T)

cmd <- "cut -f1,3 /home/guelfi/eQTL/465lymphoblastoid/EUR373.gene.cis.FDR5.best.rs137.txt | cut -d '.' -f1 | cut -f1,2"
lym.eQTLs.best.gene <- read.delim(pipe(cmd),header=F)
lym.eQTLs.best <- read.delim(file="/home/guelfi/eQTL/465lymphoblastoid/EUR373.gene.cis.FDR5.best.rs137.txt",header=F)

length(intersect(unique(eQTLPUTM$gene),unique(lym.eQTLs.gene$GENE_ID)))

length(intersect(unique(eQTLPUTM$gene),unique(lym.eQTLs.best.gene[,2])))



############################
#### Comparison Networks ###
############################

load(file="data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
load(file="data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")

netwPUTM <- read.csv("data/networs/network1237918.4.PUTM.6.rds_MM.csv")
netwSNIG <- read.csv("data/networs/network1237511.4.SNIG.6.rds_MM.csv")


eQTLSNIG <- eQTLSNIG[which(eQTLSNIG$myFDR<0.01),]
eQTLPUTM <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]


percentage <- (table(netwPUTM$module)/nrow(netwPUTM))*100
barplot(sort(table(netwPUTM$module),decreasing = T),las=2,col = names(sort(table(netwPUTM$module),decreasing = T)) )
barplot(sort(percentage,decreasing = T),las=2,col = names(sort(percentage,decreasing = T)),
        main= "Percentage of genes for each module PUTM") 

#PUTM
modules <- netwPUTM[which(netwPUTM$gene %in% eQTLPUTM$gene),"module"]
barplot(sort(table(modules),decreasing = T),las=2,col = names(sort(table(modules),decreasing = T)) )
## checking for correlation between streghness in the eQTL and hub modules 
eQTLmod <- netwPUTM[which(netwPUTM$gene %in% eQTLPUTM$gene),]
cor(eQTLmod[,"quantile"],eQTLPUTM[as.character(eQTLmod$gene),"myFDR"])

percentage <- (table(eQTLmod$module)/nrow(eQTLmod))*100
barplot(sort(percentage,decreasing = T),las=2,col = names(sort(percentage,decreasing = T)),
        main= "Percentage of eQTLs for each module PUTM") 

#SNIG
percentage <- (table(netwSNIG$module)/nrow(netwSNIG))*100
barplot(sort(percentage,decreasing = T),las=2,col = names(sort(percentage,decreasing = T)),
        main= "Percentage of genes for each module SNIG") 

modules <- netwSNIG[which(netwSNIG$ensgene %in% eQTLSNIG$gene),"module"]
barplot(sort(table(modules),decreasing = T),las=2,col = names(sort(table(modules),decreasing = T)) )
## checking for correlation between streghness in the eQTL and hub modules 
eQTLmod <- netwSNIG[which(netwSNIG$ensgene %in% eQTLSNIG$gene),]
cor(eQTLmod[,"mm"],eQTLSNIG[as.character(eQTLmod$ensgene),"myFDR"])

percentage <- (table(eQTLmod$module)/nrow(eQTLmod))*100
barplot(sort(percentage,decreasing = T),las=2,col = names(sort(percentage,decreasing = T)),
        main= "Percentage of eQTLs for each module SNIG")


### check weather genes targetee by eQTL are more hubbier


eQTLmod <- netwPUTM[which(netwPUTM$gene %in% eQTLPUTM$gene),]


hubsGenes <- eQTLmod[which(eQTLmod$quantile>0.9),]
percentage <- (table(hubsGenes[,"module"])/nrow(hubsGenes))*100
barplot(sort(percentage,decreasing = T),las=2,col = names(sort(percentage,decreasing = T)),
        main= "Percentage of eQTLs for hub genes PUTM") 


eQTLmod <- netwSNIG[which(netwSNIG$ensgene %in% eQTLSNIG$gene),]

hist(netwSNIG$mm)
hubsGenes <- eQTLmod[which(eQTLmod$mm>0.9),]
percentage <- (table(hubsGenes[,"module"])/nrow(hubsGenes))*100
barplot(sort(percentage,decreasing = T),las=2,col = names(sort(percentage,decreasing = T)),
        main= "Percentage of eQTLs for hub genes SNIG") 









######################
### Plots bio Type ###
######################


load("data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
eQTLPUTMEI <- eQTLPUTM
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTMI <- eQTLPUTM
load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")
eQTLPUTMInter <- eQTLPUTM
eQTLPUTMExons <- read.delim("data/results/finaleQTLs/exons.PUTM.txt")
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTMEx <- eQTLPUTM
rm(eQTLPUTM)

par(mar=c(9,5,3,2))
barplot(sort((table(eQTLPUTMEx$gene_biotype)/nrow(eQTLPUTMEx)*100),decreasing = T),las=2,main = "Bio type eQTLs PUTM (whole-gene exonic)",
        col=1:length(unique(eQTLPUTMEx$gene_biotype)),ylim=c(0,70))
barplot(sort((table(eQTLPUTMI$gene_biotype)/nrow(eQTLPUTMI)*100),decreasing = T),las=2,main = "Bio type eQTLs PUTM (whole-gene intronic)",
        col=1:length(unique(eQTLPUTMI$gene_biotype)),ylim=c(0,70))
par(mar=c(11,5,3,2))
barplot(sort((table(eQTLPUTMEI$gene_biotype)/nrow(eQTLPUTMEI)*100),decreasing = T),las=2,main = "Bio type eQTLs PUTM (whole-gene exonic+intronic)",
        col=1:length(unique(eQTLPUTMEI$gene_biotype)),ylim=c(0,60))


### we check the common between protein coding and intronic

common <- intersect(eQTLPUTMEx$gene,eQTLPUTMI$gene) 
tmp <- intersect(eQTLPUTMEI$gene,eQTLPUTMI$gene)
## get the specific eQTLs for intronic and exonic
intronicExSpe <- common[which(!common%in%tmp)]
 
eQTLPUTMEx[which(eQTLPUTMEx$gene %in% intronicExSpe),]
eQTLPUTMI[which(eQTLPUTMI$gene %in% intronicExSpe),]


gene <- intronicExSpe[9]
eQTLPUTMEx[gene,]; eQTLPUTMI[gene,] 

######################
### intersect SNps ###
######################


length(intersect(eQTLPUTMEI$gene,eQTLPUTMI$gene))

library(gplots)
venn(list(PUTMExonInt=unique(eQTLPUTMEI$gene),
          PUTMIntr=unique(eQTLPUTMI$gene)))

venn(list(PUTMExonInt=paste0(eQTLPUTMEI$snps,eQTLPUTMEI$gene),
          PUTMIntr=paste0(eQTLPUTMI$snps,eQTLPUTMI$gene)))

venn(list(PUTMExon=paste0(eQTLPUTMEx$snps,eQTLPUTMEx$gene),
          PUTMIntr=paste0(eQTLPUTMI$snps,eQTLPUTMI$gene)))

com <- intersect(paste0(eQTLPUTMEx$snps,eQTLPUTMEx$gene),
paste0(eQTLPUTMI$snps,eQTLPUTMI$gene))


head(com)



##########################################################
### plot the eQTLs by tissue accross different tissues ###
##########################################################


load("data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
eQTLPUTMEI <- eQTLPUTM
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTMI <- eQTLPUTM
load("data/results/finaleQTLs/intergenic.Ann.PUTM.rda")
eQTLPUTMInter <- eQTLPUTM
eQTLPUTMExons <- read.delim("data/results/finaleQTLs/exons.PUTM.txt",sep=" ")
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTMEx <- eQTLPUTM
rm(eQTLPUTM)



numeQTLs10 <- c(length(unique(eQTLPUTMEI$gene)),
                length(unique(eQTLPUTMEx$gene)),
                length(unique(eQTLPUTMI$gene)),
                length(unique(eQTLPUTMInter$gene)),
                length(unique(eQTLPUTMExons$gene)))

                
names(numeQTLs10) <- c("GeneExo+Int","GeneExonic","GeneIntronic","Intergenic","Exons")

barplot(numeQTLs10,col='red',main= "eQTLs PUTM",border=F,
        sub=paste("PUTM"))

numeQTLs5 <- c(length(unique(eQTLPUTMEI[which(eQTLPUTMEI$myFDR<0.05),"gene"])),
               length(unique(eQTLPUTMEx[which(eQTLPUTMEx$myFDR<0.05),"gene"])),
               length(unique(eQTLPUTMI[which(eQTLPUTMI$myFDR<0.05),"gene"])),
               length(unique(eQTLPUTMInter[which(eQTLPUTMInter$myFDR<0.05),"gene"])),
               length(unique(eQTLPUTMExons[which(eQTLPUTMExons$myFDR<0.05),"gene"])))

names(numeQTLs5) <- c("GeneExo+Int","GeneExonic","GeneIntronic","Intergenic","Exons")
barplot(numeQTLs5,add=T,col='orange',border=F)

numeQTLs1 <- c(length(unique(eQTLPUTMEI[which(eQTLPUTMEI$myFDR<0.01),"gene"])),
               length(unique(eQTLPUTMEx[which(eQTLPUTMEx$myFDR<0.01),"gene"])),
               length(unique(eQTLPUTMI[which(eQTLPUTMI$myFDR<0.01),"gene"])),
               length(unique(eQTLPUTMInter[which(eQTLPUTMInter$myFDR<0.01),"gene"])),
               length(unique(eQTLPUTMExons[which(eQTLPUTMExons$myFDR<0.01),"gene"])))

names(numeQTLs1) <- c("GeneExo+Int","GeneExonic","GeneIntronic","Intergenic","Exons")
barplot(numeQTLs1,add=T,col='darkgreen',border=F)
legend("topleft",c("10%","5%","1%"),col=c('red','orange','darkgreen'),title="FDR",pch=15)

## table of the independecy of the signals
sign <- cbind(table(eQTLPUTM$degree),table(eQTLSNIG$degree))

par(mar=c(9,3,1,1))
barplot(sort(table(eQTLPUTM$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLPUTM$gene_biotype))),las=2)

barplot(sort(table(eQTLSNIG$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLSNIG$gene_biotype))),las=2,main="test")

rm(numeQTLs1,numeQTLs10,numeQTLs5)



############
### SNIG ###
############


load("data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
eQTLSNIGEI <- eQTLSNIG
load("data/results/finaleQTLs/intronic.Ann.SNIG.rda")
eQTLSNIGI <- eQTLSNIG
load("data/results/finaleQTLs/intergenic.Ann.SNIG.rda")
eQTLSNIGInter <- eQTLSNIG
eQTLSNIGExons <- read.delim("data/results/finaleQTLs/exons.SNIG.txt",sep=" ")
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
eQTLSNIGEx <- eQTLSNIG
rm(eQTLSNIG)



numeQTLs10 <- c(length(unique(eQTLSNIGEI$gene)),
                length(unique(eQTLSNIGEx$gene)),
                length(unique(eQTLSNIGI$gene)),
                length(unique(eQTLSNIGInter$gene)),
                length(unique(eQTLSNIGExons$gene)))


names(numeQTLs10) <- c("GeneExo+Int","GeneExonic","GeneIntronic","Intergenic","Exons")

barplot(numeQTLs10,col='red',main= "eQTLs SNIG",border=F,
        sub=paste("SNIG"))

numeQTLs5 <- c(length(unique(eQTLSNIGEI[which(eQTLSNIGEI$myFDR<0.05),"gene"])),
               length(unique(eQTLSNIGEx[which(eQTLSNIGEx$myFDR<0.05),"gene"])),
               length(unique(eQTLSNIGI[which(eQTLSNIGI$myFDR<0.05),"gene"])),
               length(unique(eQTLSNIGInter[which(eQTLSNIGInter$myFDR<0.05),"gene"])),
               length(unique(eQTLSNIGExons[which(eQTLSNIGExons$myFDR<0.05),"gene"])))

names(numeQTLs5) <- c("GeneExo+Int","GeneExonic","GeneIntronic","Intergenic","Exons")
barplot(numeQTLs5,add=T,col='orange',border=F)

numeQTLs1 <- c(length(unique(eQTLSNIGEI[which(eQTLSNIGEI$myFDR<0.01),"gene"])),
               length(unique(eQTLSNIGEx[which(eQTLSNIGEx$myFDR<0.01),"gene"])),
               length(unique(eQTLSNIGI[which(eQTLSNIGI$myFDR<0.01),"gene"])),
               length(unique(eQTLSNIGInter[which(eQTLSNIGInter$myFDR<0.01),"gene"])),
               length(unique(eQTLSNIGExons[which(eQTLSNIGExons$myFDR<0.01),"gene"])))

names(numeQTLs1) <- c("GeneExo+Int","GeneExonic","GeneIntronic","Intergenic","Exons")
barplot(numeQTLs1,add=T,col='darkgreen',border=F)
legend("topleft",c("10%","5%","1%"),col=c('red','orange','darkgreen'),title="FDR",pch=15)

## table of the independecy of the signals
sign <- cbind(table(eQTLSNIG$degree),table(eQTLSNIG$degree))

par(mar=c(9,3,1,1))
barplot(sort(table(eQTLSNIG$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLSNIG$gene_biotype))),las=2)

barplot(sort(table(eQTLSNIG$gene_biotype),decreasing=T),
        col=c(1:length(table(eQTLSNIG$gene_biotype))),las=2,main="test")

rm(numeQTLs1,numeQTLs10,numeQTLs5)





########################
### Mina Suggestion ####
########################

## PUTM

load("data/results/finaleQTLs/exonicIntronic.Ann.SNIG.rda")
eQTLSNIGEI <- eQTLSNIG
load("data/results/finaleQTLs/intronic.Ann.SNIG.rda")
eQTLSNIGI <- eQTLSNIG
load("data/results/finaleQTLs/geneExonic.Ann.SNIG.rda")
eQTLSNIGEx <- eQTLSNIG
rm(eQTLSNIG)

library(gplots)
## overlap exonic and intronic

venn(list(SNIGExon=paste0(eQTLSNIGEx$snps,eQTLSNIGEx$external_gene_id),
          SNIGIntr=paste0(eQTLSNIGI$snps,eQTLSNIGI$external_gene_id)))

## in common we have 127 eQTLs

commonEQTLs <- intersect(paste0(eQTLSNIGEx$snps,";",eQTLSNIGEx$external_gene_id),
                    paste0(eQTLSNIGI$snps,";",eQTLSNIGI$external_gene_id))

## we get the exonic specific eQTLs
ExSpe <- paste0(eQTLSNIGEx$snps,";",eQTLSNIGEx$external_gene_id)
ExSpe <- ExSpe[-which(ExSpe %in% commonEQTLs)]
library (plyr)

geneSpe <- unlist(lapply(strsplit(ExSpe,";")
                          ,function(x){x[2]}))


##install.packages("/Users/mguelfi/Downloads/gProfileR_0.5.3.tar.gz", repos = NULL, type="source")
require(gProfileR)


go <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLSNIGI$external_gene_id,eQTLSNIGEx$external_gene_id), src_filter = "GO")
##write.csv(unique(geneSpe),file="tmp/geneSpe.csv")
write.csv(go,file="data/results/GO/GOExonSpeSNIG.csv")

InSpe <- paste0(eQTLSNIGI$snps,";",eQTLSNIGI$external_gene_id)
InSpe <- InSpe[-which(InSpe %in% commonEQTLs)]
geneSpe <- unlist(lapply(strsplit(InSpe,":")
                         ,function(x){x[2]}))

go <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLSNIGI$external_gene_id,eQTLSNIGEx$external_gene_id), src_filter = "GO")

## no point to save the file since there is no significant term

## PUTM

load("data/results/finaleQTLs/exonicIntronic.Ann.PUTM.rda")
eQTLPUTMEI <- eQTLPUTM
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTMI <- eQTLPUTM
load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTMEx <- eQTLPUTM
rm(eQTLPUTM)

library(gplots)
## overlap exonic and intronic

venn(list(PUTMExon=paste0(eQTLPUTMEx$snps,eQTLPUTMEx$external_gene_id),
          PUTMIntr=paste0(eQTLPUTMI$snps,eQTLPUTMI$external_gene_id)))

## in common we have 127 eQTLs

commonEQTLs <- intersect(paste0(eQTLPUTMEx$snps,";",eQTLPUTMEx$external_gene_id),
                         paste0(eQTLPUTMI$snps,";",eQTLPUTMI$external_gene_id))

## we get the exonic specific eQTLs
ExSpe <- paste0(eQTLPUTMEx$snps,";",eQTLPUTMEx$external_gene_id)
ExSpe <- ExSpe[-which(ExSpe %in% commonEQTLs)]
library (plyr)

geneSpe <- unlist(lapply(strsplit(ExSpe,";")
                         ,function(x){x[2]}))


##install.packages("/Users/mguelfi/Downloads/gProfileR_0.5.3.tar.gz", repos = NULL, type="source")
require(gProfileR)


go <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLPUTMI$external_gene_id,eQTLPUTMEx$external_gene_id), src_filter = "GO")
##write.csv(unique(geneSpe),file="tmp/geneSpe.csv")
write.csv(go,file="data/results/GO/GOExonSpePUTM.csv")

InSpe <- paste0(eQTLPUTMI$snps,";",eQTLPUTMI$external_gene_id)
InSpe <- InSpe[-which(InSpe %in% commonEQTLs)]
geneSpe <- unlist(lapply(strsplit(InSpe,":")
                         ,function(x){x[2]}))

go <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLPUTMI$external_gene_id,eQTLPUTMEx$external_gene_id), src_filter = "GO")
head(go)




