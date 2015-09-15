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









