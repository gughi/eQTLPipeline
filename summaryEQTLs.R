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

eQTLSNIGEx <- eQTLSNIGEx[which(eQTLSNIGEx$myFDR<0.01),]
eQTLSNIGI <- eQTLSNIGI[which(eQTLSNIGI$myFDR<0.01),]

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


go <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLSNIGI$external_gene_id,eQTLSNIGEx$external_gene_id))
##write.csv(unique(geneSpe),file="tmp/geneSpe.csv")
write.csv(go,file="data/results/GO/GOExonSpeSNIG.csv")

InSpe <- paste0(eQTLSNIGI$snps,";",eQTLSNIGI$external_gene_id)
InSpe <- InSpe[-which(InSpe %in% commonEQTLs)]
geneSpe <- unlist(lapply(strsplit(InSpe,";")
                         ,function(x){x[2]}))

go <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLSNIGI$external_gene_id,eQTLSNIGEx$external_gene_id))
write.csv(go,file="data/results/GO/GOIntronSpeSNIG.csv")

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
geneSpe <- unlist(lapply(strsplit(InSpe,";")
                         ,function(x){x[2]}))

go <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLPUTMI$external_gene_id,eQTLPUTMEx$external_gene_id), src_filter = "GO")
head(go)


kegg <- gprofiler(geneSpe ,correction_method="bonferroni",custom_bg = c(eQTLPUTMI$external_gene_id,eQTLPUTMEx$external_gene_id), src_filter = "KEGG")
head(kegg)

library("biomaRt")
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")

entrezID <- getBM(attributes=c("entrezgene"),
                  verbose = T,
                  filters="ensembl_gene_id",
                  values=geneSpe, mart=ensembl)

DOres <- enrichDO(entrezID[,1], pvalueCutoff=0.01)
head(summary(DOres))


#####################################
###  OVerlap unsentinalised eQTLs ###
#####################################



eQTLPUTMI <- read.delim(file="data/results/finaleQTLs/intronic.unsentinalised.PUTM.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLPUTMI$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
rownames(df) <- NULL
eQTLPUTMIun <- df 
colnames(eQTLPUTMIun) <- c("snps","gene","statistic","pvalue","FDR","beta","myFDR")
rm(df)

eQTLSNIGI <- read.delim(file="data/results/finaleQTLs/intronic.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLSNIGI$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
rownames(df) <- NULL
head(df)
eQTLSNIGIun <- df 
colnames(eQTLSNIGIun) <- c("snps","gene","statistic","pvalue","FDR","beta","myFDR")
save(eQTLSNIGIun,eQTLPUTMIun,file="data/results/finaleQTLs/intronic.unsentinalised.rda")
rm(df)

eQTLPUTMEx <- read.delim(file="data/results/finaleQTLs/eQTLgeneExons.unsentinalised.PUTM.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLPUTMEx$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
rownames(df) <- NULL
eQTLPUTMExun <- df 
colnames(eQTLPUTMExun) <- c("snps","gene","statistic","pvalue","FDR","beta","myFDR")
rm(df)

eQTLSNIGEx <- read.delim(file="data/results/finaleQTLs/eQTLgeneExons.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLSNIGEx$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
rownames(df) <- NULL
head(df)
eQTLSNIGExun <- df 
colnames(eQTLSNIGExun) <- c("snps","gene","statistic","pvalue","FDR","beta","myFDR")
save(eQTLSNIGExun,eQTLPUTMExun,file="data/results/finaleQTLs/geneExonic.unsentinalised.rda")
rm(df)


load("data/results/finaleQTLs/geneExonic.unsentinalised.rda")
load("data/results/finaleQTLs/intronic.unsentinalised.rda")


# filter by FDR 0.01
eQTLPUTMExun <- eQTLPUTMExun[which(as.numeric(eQTLPUTMExun[,7]) < 0.01),]
eQTLPUTMIun <- eQTLPUTMIun[which(as.numeric(eQTLPUTMIun[,7]) < 0.01),]

library(gplots)
venn(list(PUTMExon=paste0(eQTLPUTMExun[,1],eQTLPUTMExun[,2]),
          PUTMIntr=paste0(eQTLPUTMIun[,1],eQTLPUTMIun[,2])))




commonEQTLs <- intersect(paste0(eQTLPUTMExun[,1],";",eQTLPUTMExun[,2]),
                         paste0(eQTLPUTMIun[,1],";",eQTLPUTMIun[,2]))


commonGenes <- unlist(lapply(strsplit(commonEQTLs,";")
                         ,function(x){x[2]}))


## only for the plot
venn(list(PUTMExon=1:(length(unique(eQTLPUTMExun[,2]))),
          PUTMIntr=(length(unique(eQTLPUTMExun[,2]))-(length(unique(commonGenes)))+1):(length(unique(eQTLPUTMExun[,2]))+(length(unique(eQTLPUTMIun[,2]))-length(unique(commonGenes))))))

ExSpe <- unique(eQTLPUTMExun[,2])[-which(unique(eQTLPUTMExun[,2]) %in% commonGenes)]
require(gProfileR)

go <- gprofiler(ExSpe ,correction_method="bonferroni",
                custom_bg = c(unique(eQTLPUTMExun[,2],unique(eQTLPUTMIun[,2]))))

write.csv(go,file="data/results/GO/GOExonSpePUTM.csv")


InSpe <- unique(eQTLPUTMIun[,2])[-which(unique(eQTLPUTMIun[,2]) %in% commonGenes)]

go <- gprofiler(InSpe ,correction_method="bonferroni",
                custom_bg = c(unique(eQTLPUTMExun[,2],unique(eQTLPUTMIun[,2]))))

write.csv(go,file="data/results/GO/GOIntSpePUTM.csv")

###### SNIG #######

eQTLSNIGExun <- eQTLSNIGExun[which(as.numeric(eQTLSNIGExun[,7]) < 0.01),]
eQTLSNIGIun <- eQTLSNIGIun[which(as.numeric(eQTLSNIGIun[,7]) < 0.01),]


venn(list(SNIGExon=paste0(eQTLSNIGExun[,1],eQTLSNIGExun[,2]),
          SNIGIntr=paste0(eQTLSNIGIun[,1],eQTLSNIGIun[,2])))




commonEQTLs <- intersect(paste0(eQTLSNIGExun[,1],";",eQTLSNIGExun[,2]),
                         paste0(eQTLSNIGIun[,1],";",eQTLSNIGIun[,2]))


commonGenes <- unlist(lapply(strsplit(commonEQTLs,";")
                             ,function(x){x[2]}))


## only for the plot
venn(list(SNIGExon=1:(length(unique(eQTLSNIGExun[,2]))),
          SNIGIntr=(length(unique(eQTLSNIGExun[,2]))-(length(unique(commonGenes)))+1):(length(unique(eQTLSNIGExun[,2]))+(length(unique(eQTLSNIGIun[,2]))-length(unique(commonGenes))))))


ExSpe <- unique(eQTLSNIGExun[,2])[-which(unique(eQTLSNIGExun[,2]) %in% commonGenes)]
require(gProfileR)

go <- gprofiler(ExSpe ,correction_method="bonferroni",
                custom_bg = c(unique(eQTLSNIGExun[,2],unique(eQTLSNIGIun[,2]))))

write.csv(go,file="data/results/GO/GOExonSpeSNIG.csv")


InSpe <- unique(eQTLSNIGIun[,2])[-which(unique(eQTLSNIGIun[,2]) %in% commonGenes)]

go <- gprofiler(InSpe ,correction_method="bonferroni",
                custom_bg = c(unique(eQTLSNIGExun[,2],unique(eQTLSNIGIun[,2]))))

write.csv(go,file="data/results/GO/GOIntSpeSNIG.csv")

####################################################
#### Overlap exonic,intronic and exonic-intronic ###
####################################################

eQTLPUTMEI <- read.delim(file="data/results/finaleQTLs/eQTLExonIntrons.unsentinalised.PUTM.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLPUTMEI$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
rownames(df) <- NULL
eQTLPUTMEIun <- df 
colnames(eQTLPUTMEIun) <- c("snps","gene","statistic","pvalue","FDR","beta","myFDR")
rm(df)

eQTLSNIGEI <- read.delim(file="data/results/finaleQTLs/eQTLExonIntrons.unsentinalised.SNIG.txt",  as.is=T, header=T)
tmp <- strsplit(eQTLSNIGEI$SNP," ")
df <- ldply(list(tmp), data.frame)
rm(tmp)
df <- t(df)
rownames(df) <- NULL
head(df)
eQTLSNIGEIun <- df 
colnames(eQTLSNIGEIun) <- c("snps","gene","statistic","pvalue","FDR","beta","myFDR")
save(eQTLSNIGEIun,eQTLPUTMEIun,file="data/results/finaleQTLs/eQTLExonIntrons.unsentinalised.rda")
rm(df)

load("data/results/finaleQTLs/geneExonic.unsentinalised.rda")
load("data/results/finaleQTLs/intronic.unsentinalised.rda")
load("data/results/finaleQTLs/eQTLExonIntrons.unsentinalised.rda")

eQTLPUTMExun <- eQTLPUTMExun[which(as.numeric(eQTLPUTMExun[,7]) < 0.01),]
eQTLPUTMIun <- eQTLPUTMIun[which(as.numeric(eQTLPUTMIun[,7]) < 0.01),]
eQTLPUTMEIxun <- eQTLPUTMExun[which(as.numeric(eQTLPUTMExun[,7]) < 0.01),]

## common eQTL for the three groups
commonEQTLs <- intersect(intersect(paste0(eQTLSNIGExun[,1],";",eQTLSNIGExun[,2]),
                         paste0(eQTLSNIGIun[,1],";",eQTLSNIGIun[,2])),
                         paste0(eQTLSNIGEIun[,1],";",eQTLSNIGEIun[,2]))


commonGenes <- unlist(lapply(strsplit(commonEQTLs,";")
                             ,function(x){x[2]}))

## intersection between the three groups
length(unique(commonGenes))
# 171

commonEQTLs <- intersect(paste0(eQTLSNIGExun[,1],";",eQTLSNIGExun[,2]),
                                   paste0(eQTLSNIGIun[,1],";",eQTLSNIGIun[,2]))


commonGenes <- unlist(lapply(strsplit(commonEQTLs,";")
                             ,function(x){x[2]}))

## intersection between the exon -introns groups
length(unique(commonGenes))
# 173

commonEQTLs <- intersect(paste0(eQTLSNIGExun[,1],";",eQTLSNIGExun[,2]),
                         paste0(eQTLSNIGEIun[,1],";",eQTLSNIGEIun[,2]))


commonGenes <- unlist(lapply(strsplit(commonEQTLs,";")
                             ,function(x){x[2]}))

## intersection between the Exon - Exon+introns groups
length(unique(commonGenes))
# 321

commonEQTLs <- intersect(paste0(eQTLSNIGEIun[,1],";",eQTLSNIGEIun[,2]),
                         paste0(eQTLSNIGIun[,1],";",eQTLSNIGIun[,2]))


commonGenes <- unlist(lapply(strsplit(commonEQTLs,";")
                             ,function(x){x[2]}))

## intersection between the introns - Exon+introns groups
length(unique(commonGenes))
# 602



length(unique(eQTLSNIGEIun[,2]))
length(unique(eQTLSNIGIun[,2]))
length(unique(eQTLSNIGExun[,2]))




########################
## Annotation of eQTL ##
########################

load("data/results/finaleQTLs/geneExonic.unsentinalised.rda")

library(doParallel)
library(foreach)

library(biomaRt)
ensembl <- useMart(biomart="ENSEMBL_MART_SNP", host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_snp")


detectCores()
## [1] 24
# create the cluster with the functions needed to run
cl <- makeCluster(20)
clusterExport(cl, c("annSinSNP","getBM"))

registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
conse <- foreach(i=1:nrow(eQTLPUTMExun),.combine=rbind,.verbose=F)%dopar%annSinSNP(eQTLPUTMExun[i,1],ensembl)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)
colnames(conse) <- c("rsID","consequence")
save(conse,eQTLPUTMExun,file="data/results/finaleQTLs/geneExonic.un.SNPAnn.PUTM.rda")


load("data/results/finaleQTLs/intronic.unsentinalised.rda")

cl <- makeCluster(20)
clusterExport(cl, c("annSinSNP","getBM"))

registerDoParallel(cl)
getDoParWorkers()

start <- Sys.time()
conse <- foreach(i=1:nrow(eQTLPUTMIun),.combine=rbind,.verbose=F)%dopar%annSinSNP(eQTLPUTMIun[i,1],ensembl)
##exonicRegions <- foreach(i=1:20,.combine=rbind,.verbose=F)%dopar%getRegionsBED(geneIDs[i],exonsdef)
end <- Sys.time()
end-start
stopCluster(cl)
rm(cl,end,start)
colnames(conse) <- c("rsID","consequence")
save(conse,eQTLPUTMIun,file="data/results/finaleQTLs/intronic.un.SNPAnn.PUTM.rda")

## overlap of the unsentinalised eQTLs based on teh SNP annotation



load("data/results/finaleQTLs/intronic.un.SNPAnn.PUTM.rda")
# filter by FDR 0.01
eQTLPUTMIun <- cbind(eQTLPUTMIun,conse)
eQTLPUTMIun <- eQTLPUTMIun[which(as.numeric(eQTLPUTMIun[,7]) < 0.01),]
rownames(eQTLPUTMIun) <- 1:nrow(eQTLPUTMIun)
rm(conse)

load("data/results/finaleQTLs/geneExonic.un.SNPAnn.PUTM.rda")
# filter by FDR 0.01
eQTLPUTMExun <- cbind(eQTLPUTMExun,conse)
eQTLPUTMExun <- eQTLPUTMExun[which(as.numeric(eQTLPUTMExun[,7]) < 0.01),]
rownames(eQTLPUTMExun) <- 1:nrow(eQTLPUTMExun)
rm(conse)

common <-intersect(paste0(eQTLPUTMExun[,1],eQTLPUTMExun[,2]),
          paste0(eQTLPUTMIun[,1],eQTLPUTMIun[,2]))

## we plot a barplot for teh SNP annotation for teh exon Specific 
speEx <- eQTLPUTMExun[-which(paste0(eQTLPUTMExun[,1],eQTLPUTMExun[,2]) %in% common),9]
speEx <- sapply(speEx,function(x){unlist(strsplit(x,";"))[1]})
nOfNA <- length(which(speEx=="NA"))
speEx <-speEx[-which(speEx=="NA")]
par(mar=c(11,5,3,2))
barplot(sort((table(speEx)/length(speEx))*100,decreasing = T),las=2,main=paste("Common SNPs annotated gene-exons (SNPs with non annotated",nOfNA,")" ),
        col=1:length(table(speEx)),ylim=c(0,60),ylab="percentage")


## we plot a barplot for teh SNP annotation for the introns specific
speIn <- eQTLPUTMIun[-which(paste0(eQTLPUTMIun[,1],eQTLPUTMIun[,2]) %in% common),9]
speIn <- sapply(speIn,function(x){unlist(strsplit(x,";"))[1]})
nOfNA <- length(which(speIn=="NA"))
speIn <-speIn[-which(speIn=="NA")]
par(mar=c(11,5,3,2))
barplot(sort((table(speIn)/length(speIn))*100,decreasing = T),las=2,main=paste("Common SNPs annotated gene-introns (SNPs with non annotated",nOfNA,")" ),
        col=1:length(table(speIn)),ylim=c(0,60),ylab="percentage")

commonAnn <- eQTLPUTMExun[which(paste0(eQTLPUTMExun[,1],eQTLPUTMExun[,2]) %in% common),9]
commonAnn <- sapply(commonAnn,function(x){unlist(strsplit(x,";"))[1]})

## remove NAs
nOfNA <- length(which(commonAnn=="NA"))
commonAnn <-commonAnn[-which(commonAnn=="NA")]
par(mar=c(11,5,3,2))
barplot(sort((table(commonAnn)/length(commonAnn))*100,decreasing = T),las=2,main=paste("Common SNPs annotated common (SNPs with non annotated",nOfNA,")" ),
        col=1:length(table(commonAnn)),ylim=c(0,60),ylab="percentage")

# now we test the different features
par(mfrow=c(1,2))
feature <- "intron_variant"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage")


feature <- "downstream_gene_variant"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(tmp,las=2,main=paste(feature ),
        col=1:3,ylab="percentage")


feature <- "upstream_gene_variant"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(tmp,las=2,main=paste(feature ),
        col=1:3,ylab="percentage")


feature <- "nc_transcript_variant"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage")


feature <- "NMD_transcript_variant"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage")


feature <- "feature_elongation"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage")




feature <- "3_prime_UTR_variant"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage")

feature <- "feature_truncation"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage")

feature <- "missense_variant"
common <- (table(commonAnn)[feature]/length(commonAnn))*100
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(common,exonic,intronic)
names(tmp) <- c("common","exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage")









library(gplots)
### check the overlap of the genes between introns and exonic

load("data/expr/normalisedCounts/genic/geneExons/RPKM.cqn.PUTM")
geneExonic <- rownames(RPKM.cqn)
load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
geneIntronic <- rownames(RPKM.cqn)

venn(list(PUTMExon=geneIntronic,
          PUTMIntr=geneExonic))


5254/(5254+13758+4465)
##22.4%
13758/(5254+13758+4465)
## 58.6%
4465/(5254+13758+4465)
## 19%

load("data/results/finaleQTLs/geneExonic.unsentinalised.rda")
load("data/results/finaleQTLs/intronic.unsentinalised.rda")

eQTLPUTMExun <- eQTLPUTMExun[which(as.numeric(eQTLPUTMExun[,7]) < 0.01),]
eQTLPUTMIun <- eQTLPUTMIun[which(as.numeric(eQTLPUTMIun[,7]) < 0.01),]

commonGenes <- intersect(paste0(eQTLPUTMExun[,2]),
                         paste0(eQTLPUTMIun[,2]))

venn(list(PUTMExon=paste0(unique(eQTLPUTMExun[,2])),
          PUTMIntr=paste0(unique(eQTLPUTMIun[,2]))))


271/(271+423+469)
## 23%
469/(271+423+469)
## 40%
423/(271+423+469)
## 36%

ExSpe <- unique(eQTLPUTMExun[,2])[- which(unique(eQTLPUTMExun[,2]) %in% commonGenes)]

length(intersect(ExSpe,geneIntronic))
## 263

InSpe <- unique(eQTLPUTMIun[,2])[- which(unique(eQTLPUTMIun[,2]) %in% commonGenes)]

length(intersect(InSpe,geneExonic))
# 320 

271/(271+320+263)
# 32%
320/(271+320+263)
# 37%
263/(271+320+263)
# 31%

intersect(InSpe,geneExonic)

load("data/expr/normalisedCounts/genic/geneExons/RPKM.cqn.PUTM")
exprExonic <- RPKM.cqn[intersect(InSpe,geneExonic),] 
load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
exprIntronic <- RPKM.cqn[intersect(InSpe,geneExonic),] 

d<-density(exprIntronic)
plot(d, xlim=c(-10,10), main="Expression genes eQTL intronic only")
polygon(d, col='skyblue') 
polygon(density(exprExonic), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)

load("data/expr/normalisedCounts/genic/geneExons/RPKM.cqn.PUTM")
exprExonic <- RPKM.cqn[intersect(ExSpe,geneIntronic),] 
load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
exprIntronic <- RPKM.cqn[intersect(ExSpe,geneIntronic),] 
rm(RPKM.cqn)

d<-density(exprIntronic)
plot(d, xlim=c(-10,10), main="Expression genes eQTL exonic only")
polygon(d, col='skyblue') 
polygon(density(exprExonic), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)


load("data/expr/normalisedCounts/genic/geneExons/RPKM.cqn.PUTM")
exprExonic <- RPKM.cqn
load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
exprIntronic <- RPKM.cqn

d<-density(exprIntronic)
plot(d, xlim=c(-10,10), main="Expression All genes")
polygon(d, col='skyblue') 
polygon(density(exprExonic), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)

int <- intersect(rownames(exprExonic),rownames(exprIntronic))
exprExonic <- exprExonic[as.character(int),]
exprIntronic <- exprIntronic[as.character(int),]

d<-density(exprIntronic)
plot(d, xlim=c(-10,10), main="Expression All common genes")
polygon(d, col='skyblue') 
polygon(density(exprExonic), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)

d<-density(var(exprIntronic))
plot(d, xlim=c(0,4), main="Variance of the expression All common genes")
polygon(d, col='skyblue') 
polygon(density(var(exprExonic)), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)

commonGenes <- intersect(paste0(eQTLPUTMExun[,2]),
                         paste0(eQTLPUTMIun[,2]))

exprExonic <- exprExonic[as.character(commonGenes),]
exprIntronic <- exprIntronic[as.character(commonGenes),]

d<-density(exprIntronic)
plot(d, xlim=c(-10,10), main="Expression genes eQTL common")
polygon(d, col='skyblue') 
polygon(density(exprExonic), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)

d<-density(var(exprIntronic))
plot(d, xlim=c(0,3),  main="Variance eQTL common")
polygon(d, col='skyblue') 
polygon(density(var(exprExonic)), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)

############################################################
### Check differences in the TSS fro intronic and exonic ###
############################################################



load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTMEx <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTMI <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]

commonGenes <- intersect(paste0(eQTLPUTMEx[,2]),
                         paste0(eQTLPUTMI[,2]))

eQTLPUTMEx <- eQTLPUTMEx[-which(eQTLPUTMEx$gene %in% commonGenes),]

posSNP <- unlist(lapply(strsplit(as.character(eQTLPUTMEx$snps),":"),function(x){x[2]}))
DisGeneStart <- eQTLPUTMEx$TSS - as.integer(posSNP)
eQTLPUTMEx <- cbind(eQTLPUTMEx,DisGeneStart) 


eQTLPUTMI <- eQTLPUTMI[-which(eQTLPUTMI$gene %in% commonGenes),]

posSNP <- unlist(lapply(strsplit(as.character(eQTLPUTMI$snps),":"),function(x){x[2]}))
DisGeneStart <- eQTLPUTMI$TSS - as.integer(posSNP)
eQTLPUTMI <- cbind(eQTLPUTMI,DisGeneStart) 

par(mfrow=c(1,1))
hist(eQTLPUTMI$DisGeneStart,col='skyblue',border=F,main= "TSS Intronic vs Exonic",
     sub=paste("Intronic:",length(eQTLPUTMI$DisGeneStart),"Exonic:",length(eQTLPUTMEx$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(eQTLPUTMI$DisGeneStart,eQTLPUTMEx$DisGeneStart)$p.value),freq=FALSE,breaks = 40, ylim=c(0,7e-06))
hist(eQTLPUTMEx$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(eQTLPUTMI$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(eQTLPUTMEx$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("Intronic","Exonic"),col=c('skyblue','red'),pch=15)



## plots for the SNP annotation
{


load("data/expr/normalisedCounts/genic/geneExons/RPKM.cqn.PUTM")
exprExonic <- RPKM.cqn
load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
exprIntronic <- RPKM.cqn
rm(RPKM.cqn,PUTM,covs)
int <- intersect(rownames(exprExonic),rownames(exprIntronic))


load("data/results/finaleQTLs/geneExonic.Ann.PUTM.rda")
eQTLPUTMEx <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]
load("data/results/finaleQTLs/intronic.Ann.PUTM.rda")
eQTLPUTMI <- eQTLPUTM[which(eQTLPUTM$myFDR<0.01),]

commonGenes <- intersect(paste0(eQTLPUTMEx[,2]),
                         paste0(eQTLPUTMI[,2]))

eQTLPUTMEx <- eQTLPUTMEx[-which(eQTLPUTMEx$gene %in% commonGenes),]
eQTLPUTMEx <- eQTLPUTMEx[which(eQTLPUTMEx$gene %in% int),]


eQTLPUTMI <- eQTLPUTMI[-which(eQTLPUTMI$gene %in% commonGenes),]
eQTLPUTMI <- eQTLPUTMI[which(eQTLPUTMI$gene %in% int),]



posSNP <- unlist(lapply(strsplit(as.character(eQTLPUTMEx$snps),":"),function(x){x[2]}))
DisGeneStart <- eQTLPUTMEx$TSS - as.integer(posSNP)
eQTLPUTMEx <- cbind(eQTLPUTMEx,DisGeneStart) 


posSNP <- unlist(lapply(strsplit(as.character(eQTLPUTMI$snps),":"),function(x){x[2]}))
DisGeneStart <- eQTLPUTMI$TSS - as.integer(posSNP)
eQTLPUTMI <- cbind(eQTLPUTMI,DisGeneStart) 

par(mfrow=c(1,1))
hist(eQTLPUTMI$DisGeneStart,col='skyblue',border=F,main= "TSS Intronic vs Exonic",
     sub=paste("Intronic:",length(eQTLPUTMEx$DisGeneStart),"Exonic:",length(eQTLPUTMI$DisGeneStart)),
     xlab=paste("KS pvalue:",ks.test(eQTLPUTMI$DisGeneStart,eQTLPUTMEx$DisGeneStart)$p.value),freq=FALSE,breaks = 40, ylim=c(0,7e-06))
hist(eQTLPUTMEx$DisGeneStart,add=T,col=scales::alpha('red',.5),border=F,freq=FALSE,breaks=40)
lines(density(eQTLPUTMI$DisGeneStart, adjust = 2), col = "skyblue")
lines(density(eQTLPUTMEx$DisGeneStart, adjust = 2), col = "red")
legend("topright",c("Intronic","Exonic"),col=c('skyblue','red'),pch=15)



load("data/results/finaleQTLs/intronic.un.SNPAnn.PUTM.rda")
# filter by FDR 0.01
eQTLPUTMIun <- cbind(eQTLPUTMIun,conse)
eQTLPUTMIun <- eQTLPUTMIun[which(as.numeric(eQTLPUTMIun[,7]) < 0.01),]
rownames(eQTLPUTMIun) <- 1:nrow(eQTLPUTMIun)
rm(conse)

load("data/results/finaleQTLs/geneExonic.un.SNPAnn.PUTM.rda")
# filter by FDR 0.01
eQTLPUTMExun <- cbind(eQTLPUTMExun,conse)
eQTLPUTMExun <- eQTLPUTMExun[which(as.numeric(eQTLPUTMExun[,7]) < 0.01),]
rownames(eQTLPUTMExun) <- 1:nrow(eQTLPUTMExun)
rm(conse)

common <-intersect(paste0(eQTLPUTMExun[,1],eQTLPUTMExun[,2]),
                   paste0(eQTLPUTMIun[,1],eQTLPUTMIun[,2]))

## we plot a barplot for teh SNP annotation for teh exon Specific 
speEx <- eQTLPUTMExun[which(eQTLPUTMExun[,2] %in% commonGenes),]
speEx <- speEx[-which(paste0(speEx[,1],speEx[,2]) %in% common),9]

speEx <- sapply(speEx,function(x){unlist(strsplit(x,";"))[1]})
nOfNA <- length(which(speEx=="NA"))
speEx <-speEx[-which(speEx=="NA")]
par(mar=c(11,5,3,2))
barplot(sort((table(speEx)/length(speEx))*100,decreasing = T),las=2,main=paste("Common SNPs annotated gene-exons (SNPs with non annotated",nOfNA,")" ),
        col=1:length(table(speEx)),ylim=c(0,60),ylab="percentage")


## we plot a barplot for teh SNP annotation for the introns specific
speIn<- eQTLPUTMIun[which(eQTLPUTMIun[,2] %in% commonGenes),]
speIn <- speIn[-which(paste0(speIn[,1],speIn[,2]) %in% common),9]
speIn <- sapply(speIn,function(x){unlist(strsplit(x,";"))[1]})
nOfNA <- length(which(speIn=="NA"))
speIn <-speIn[-which(speIn=="NA")]
par(mar=c(11,5,3,2))
barplot(sort((table(speIn)/length(speIn))*100,decreasing = T),las=2,main=paste("Common SNPs annotated gene-introns (SNPs with non annotated",nOfNA,")" ),
        col=1:length(table(speIn)),ylim=c(0,60),ylab="percentage")


par(mfrow=c(1,3))
feature <- "intron_variant"
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(exonic,intronic)
names(tmp) <- c("exonic","intronic")
barplot(tmp,las=2,main=paste(feature ),
        col=1:3,ylab="percentage",ylim = c(0,60),sub = paste("chi-square test:",
                                                             chisq.test(c(table(speEx)[feature],table(speIn)[feature]))$p.value))


feature <- "downstream_gene_variant"
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(exonic,intronic)
names(tmp) <- c("exonic","intronic")
barplot(tmp,las=2,main=paste(feature ),
        col=1:3,ylab="percentage",ylim = c(0,16),sub = paste("chi-square test:",
                                                             chisq.test(c(table(speEx)[feature],table(speIn)[feature]))$p.value))


feature <- "upstream_gene_variant"
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(exonic,intronic)
names(tmp) <- c("exonic","intronic")
barplot(tmp,las=2,main=paste(feature ),
        col=1:3,ylab="percentage",ylim = c(0,16),sub = paste("chi-square test:",
                                                             chisq.test(c(table(speEx)[feature],table(speIn)[feature]))$p.value))


feature <- "3_prime_UTR_variant"
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(exonic,intronic)
names(tmp) <- c("exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,ylab="percentage",sub = paste("chi-square test:",
                                              chisq.test(c(table(speEx)[feature],table(speIn)[feature]))$p.value))


feature <- "3_prime_UTR_variant"
exonic <- table(speEx)[feature]
intronic <- table(speIn)[feature]
tmp <- c(exonic,intronic)
names(tmp) <- c("exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,sub = paste("Total Exonic=",length(speEx),"Intronic",length(speIn)))



par(mfrow=c(1,3))
feature <- "NMD_transcript_variant"
exonic <- (table(speEx)[feature]/length(speEx))*100
intronic <- (table(speIn)[feature]/length(speIn))*100
tmp <- c(exonic,intronic)
names(tmp) <- c("exonic","intronic")
barplot(tmp,las=2,main=paste(feature ),
        col=1:3,ylab="percentage",,sub = paste("chi-square test:",
                                                             chisq.test(c(table(speEx)[feature],table(speIn)[feature]))$p.value))

feature <- "NMD_transcript_variant"
exonic <- table(speEx)[feature]
intronic <- table(speIn)[feature]
tmp <- c(exonic,intronic)
names(tmp) <- c("exonic","intronic")
barplot(sort(tmp,decreasing = T),las=2,main=paste(feature ),
        col=1:3,sub = paste("Total Exonic=",length(speEx),"Intronic",length(speIn)))

}


load("data/expr/normalisedCounts/genic/geneExons/RPKM.cqn.PUTM")
exprExonic <- RPKM.cqn
load("data/expr/normalisedCounts/genic/geneIntronic/RPKM.cqn.PUTM")
exprIntronic <- RPKM.cqn
rm(RPKM.cqn,PUTM,covs)
int <- intersect(rownames(exprExonic),rownames(exprIntronic))

head(rownames(exprExonic)[-which(rownames(exprExonic) %in% int)])
## "ENSG00000000460" "ENSG00000000938" "ENSG00000000971" "ENSG00000001626" "ENSG00000002822" "ENSG00000003137"
head(rownames(exprIntronic)[-which(rownames(exprIntronic) %in% int)])
## "ENSG00000001630" "ENSG00000006555" "ENSG00000008516" "ENSG00000008517" "ENSG00000010704" "ENSG00000011052"



ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="Jun2013.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")
head(listAttributes(ensembl))

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","chromosome_name","start_position","end_position","gene_biotype"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=rownames(exprIntronic)[-which(rownames(exprIntronic) %in% int)], mart=ensembl)

##ENSG00000006283    CACNA1G 17  48638429  48704835 protein_coding

par(mar=c(11,5,3,2))
barplot(sort(table(geneNames$gene_biotype),decreasing = T),las=2)

geneNames <- getBM(attributes=c("ensembl_gene_id","external_gene_id","chromosome_name","start_position","end_position","gene_biotype"),
                   verbose = T,
                   filters="ensembl_gene_id",
                   values=rownames(exprExonic)[-which(rownames(exprExonic) %in% int)], mart=ensembl)



par(mar=c(11,5,3,2))
barplot(sort(table(geneNames$gene_biotype),decreasing = T),las=2)

head(geneNames,20)

##ENSG00000006555            TTC22               1       55245385     55266940 protein_coding



load_all()

## plot read-depth
plotReadDepth("ENSG00000006555")

## boxplot raw counts
load("data/expr/rawCounts/genic/exprIntrons.rda")
load("data/expr/rawCounts/genic/exprSQ.rda")
load("data/general/sampleInfo.rda")

PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
exprSQ <- exprSQ["ENSG00000006555",as.character(PUTM$A.CEL_file)]
intronicReads <- intronicReads[as.character(PUTM$A.CEL_file),"ENSG00000006555"]

par(mar=c(3,3,3,1))
boxplot(as.numeric(exprSQ[1,]),intronicReads,names=c("exonic","intronic"),main="Raw counts ENSG00000006555")


## boxplot RPKM1
load("data/general/genesWidthExonic.rda")

load("data/expr/rawCounts/genic/exprSQ.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")

## convert the genes that have NAs
exprSQ[is.na(exprSQ)]=0
## remove genes that not expressed in any gene
exprSQ <- exprSQ[rowSums(exprSQ>0)>0,]
cat("Processing PUTM \n")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]

# now we select the expression for the PUTM only samples
expr <- exprSQ[,as.character(PUTM$A.CEL_file)]
rm(exprSQ)

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[as.character(rownames(expr))]
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

library(easyRNASeq)
RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

RPKM.std <- RPKM.std["ENSG00000006555",as.character(PUTM$A.CEL_file)]


load("data/general/genesWidthExonic.rda")

load("data/expr/rawCounts/genic/exprIntrons.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")


intronicReads <- t(intronicReads)
## convert the genes that have NAs
intronicReads[is.na(intronicReads)]=0
## remove genes that not expressed in any gene
intronicReads <- intronicReads[rowSums(intronicReads>0)>0,]
cat("Processing PUTM \n")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]

# now we select the expression for the PUTM only samples
expr <- intronicReads[,as.character(PUTM$A.CEL_file)]
rm(intronicReads)

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[as.character(rownames(expr))]
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

library(easyRNASeq)
RPKM.stdIntronic <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)

RPKM.stdIntronic <- RPKM.stdIntronic["ENSG00000006555",as.character(PUTM$A.CEL_file)]


par(mar=c(3,3,3,1))
boxplot(RPKM.std,RPKM.stdIntronic,names=c("exonic","intronic"),main="RPKM ENSG00000006555")



RPKM.stdIntronic=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]



load("data/expr/rawCounts/genic/exprIntrons.rda")
load("data/expr/rawCounts/genic/exprSQ.rda")
load("data/general/sampleInfo.rda")

PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified =="PUTM"),]
exprSQ <- exprSQ["ENSG00000006555",as.character(PUTM$A.CEL_file)]
intronicReads <- intronicReads[as.character(PUTM$A.CEL_file),"ENSG00000006555"]

par(mar=c(3,3,3,1))
boxplot(as.numeric(exprSQ[1,]),intronicReads,names=c("exonic","intronic"),main="Raw counts ENSG00000006555")


## boxplot RPKM1
load("data/general/genesWidthExonic.rda")

load("data/expr/rawCounts/genic/exprSQ.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")

## convert the genes that have NAs
exprSQ[is.na(exprSQ)]=0
## remove genes that not expressed in any gene
exprSQ <- exprSQ[rowSums(exprSQ>0)>0,]
cat("Processing PUTM \n")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]

# now we select the expression for the PUTM only samples
expr <- exprSQ[,as.character(PUTM$A.CEL_file)]
rm(exprSQ)

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[as.character(rownames(expr))]
stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

library(easyRNASeq)
RPKM.std <- RPKM(as.matrix(expr), NULL, 
                 lib.size=librarySize, 
                 feature.size=length)


load("data/general/geneswidthIntornic.rda")
load("data/expr/rawCounts/genic/exprIntrons.rda")
# load the sample info to get the IDs for each tissue
load("data/general/sampleInfo.rda")

intronicReads <- t(intronicReads)
## convert the genes that have NAs
intronicReads[is.na(intronicReads)]=0
table(is.na(intronicReads))
## clean negative reads
table(intronicReads[intronicReads<0])
intronicReads[intronicReads<0]=0
table(intronicReads[intronicReads<0])

## remove genes that not expressed in any gene
intronicReads <- intronicReads[rowSums(intronicReads>0)>0,]
cat("Processing PUTM \n")
PUTM <- sampleInfo[which(sampleInfo$U.Region_simplified=="PUTM"),]

length <- as.numeric(geneswidth[,2])
names(length) <-  as.character(geneswidth[,1])
length <- length[!is.na(length)]

# now we select the expression for the PUTM only samples
expr <- intronicReads[names(length),as.character(PUTM$A.CEL_file)]
rm(intronicReads)

librarySize <- read.csv(file="data/general/librarySize.csv", row.names=1)
librarySize <- librarySize[as.character(PUTM$A.CEL_file),]
names(librarySize) <- as.character(PUTM$A.CEL_file)

stopifnot(identical(colnames(expr),names(librarySize)))
stopifnot(identical(rownames(expr),names(length)))              

library(easyRNASeq)
RPKM.stdIntronic <- RPKM(as.matrix(expr), NULL, 
                         lib.size=librarySize, 
                         feature.size=length)

## remove intronic genes that have 0 as length (probably genes composed by only exons)
genesList <- rownames(RPKM.stdIntronic)
genesList <- genesList[!is.na(genesList)]
RPKM.stdIntronic <- RPKM.stdIntronic[as.character(genesList),]
rm(genesList)

int <- intersect(rownames(RPKM.stdIntronic),rownames(RPKM.std))

d<-density(log(RPKM.stdIntronic[as.character(int),]))
plot(d, main="log(RPKM) density for all the genes without filtering")
polygon(d, col='skyblue') 
polygon(density(log(RPKM.std[as.character(int),])), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)


RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
RPKM.stdIntronic=RPKM.stdIntronic[rowSums(RPKM.stdIntronic>=0.1)>(ncol(RPKM.stdIntronic)-((ncol(RPKM.stdIntronic)*20)/100)),]


int <- intersect(rownames(RPKM.stdIntronic),rownames(RPKM.std))

d<-density(log(RPKM.stdIntronic[as.character(int),]))
plot(d, main="log(RPKM) density for all the genes filtered")
polygon(d, col='skyblue') 
polygon(density(log(RPKM.std[as.character(int),])), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)


GCconExo <- read.delim("../Bioinformatics/ensemblRef/GCcontentExonic",sep=" ",row.names=1)
load("data/general/GCcontentIntronic.rda")
rownames(GCcontentByGene) <- GCcontentByGene[,1]
GCcontentByGene <- GCcontentByGene[,-1]
GCconInt <- GCcontentByGene
rm(GCcontentByGene)

GCconExo <- GCconExo[as.character(int),]
names(GCconExo) <- as.character(int)
GCconInt <- GCconInt[as.character(int)]

d<-density(as.numeric(GCconInt))
plot(d, main="GC content")
polygon(d, col='skyblue') 
polygon(density(as.numeric(GCconExo)), col=scales::alpha('red',.5)) 
legend("topright",c("expr Intronic","expr Exonic"),col=c('skyblue','red'),pch=15)
rm(d)













