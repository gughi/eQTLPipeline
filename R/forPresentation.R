load("data/expr/normalisedCounts/intergenic/resids.PUTM.rda")
load("data/expr/normalisedCounts/intergenic/RPKM.cqn.PUTM")

hist(resids,main="residuals distribution",breaks = 40,
     sub=paste("range max",max(resids),"min",min(resids),
               "\nmedian",median(resids),"mean:",mean(resids)),xlab="")

hist(RPKM.cqn,main="residuals distribution",breaks = 40,
     sub=paste("range max",max(RPKM.cqn),"min",min(RPKM.cqn),
               "\nmedian",median(RPKM.cqn),"mean:",mean(RPKM.cqn)),xlab="")


dim(RPKM.cqn)
dim(starStopReg)
dim(resids)

# Does the RPKM corrects for the region length?
correlation <- apply(RPKM.cqn,2,function(x){cor(x,starStopReg$width)})
hist(correlation)

# Does the residualas corrects for the region length?
correlation <- apply(resids,1,function(x){cor(x,starStopReg$width)})
hist(correlation)



met <- read.csv("data/general/QCmetrics.csv")

par(mar=c(7,5,3,2))
reads <- c(mean(met$A.Reads),0,0,0)
names(reads) <-c("Total Reads","","","" )
barplot(reads,col=c(139,73,35,17),las=2,main="Reads")

par(mar=c(7,5,3,2))
reads <- c(mean(met$A.Reads),mean(met$F.TotReadsNoAdap,na.rm = T),0,0)
names(reads) <-c("Total Reads","No adapters","","")
barplot(reads,col=c(139,28,35,17),las=2,main="Reads")

bamReads <- read.delim("data/general/readsBam.txt",header = F)

par(mar=c(7,5,3,2))
reads <- c(mean(met$A.Reads),mean(met$F.TotReadsNoAdap,na.rm = T),mean(bamReads[,1]),0)
names(reads) <-c("Total Reads","No adapters","Total Mapped","")
barplot(reads,col=c(139,28,34,17),las=2,main="Reads")

par(mar=c(7,5,3,2))
reads <- c(mean(met$A.Reads),mean(met$F.TotReadsNoAdap,na.rm = T),mean(bamReads[,1]),mean(met$R.Total.NorRNA))
names(reads) <-c("Total Reads","No adapters","Total Mapped","no rRNA")
barplot(reads,col=c(139,28,34,16),las=2,main="Reads")


par(mar=c(7,5,3,2))
reads <- c(mean(met$A.Reads),mean(met$F.TotReadsNoAdap,na.rm = T),mean(bamReads[,1]),mean(met$R.Total.NorRNA),mean(met$R.Exonic.Reads))
names(reads) <-c("Total Reads","No adapters","Total Mapped","no rRNA","Exonic Reads")
barplot(reads,col=c(139,28,34,16,28),las=2,main="Reads")




mean(met$F.TotReadsNoAdap,na.rm = T)/mean(met$A.Reads)
mean(bamReads[,1])/mean(met$A.Reads)
mean(met$R.Total.NorRNA)/mean(met$A.Reads)
mean(met$R.Exonic.Reads)/mean(met$A.Reads)

## collect information of bam files:
##Apollo
####finalOut <- cbind.data.frame(transcriptID)

for (i in 1:length(met$A.CEL_file))
{
  print(paste(met$A.CEL_file[i],date()))
  
  ##open the cufflink file
  if(file.exists(paste0("/home/ukbec/brainRNAseq/tophat2_June2013/Sample_",met$A.CEL_file[i],"/trim_galore/align_summary.txt")))
  {  
    cmd <- paste0("grep 'Mapped' /home/ukbec/brainRNAseq/tophat2_June2013/Sample_",met$A.CEL_file[i],"/trim_galore/align_summary.txt | cut -f17 -d' ' | awk '{ sum+=$1} END {print sum}'  >> ~/test.txt")
    system(cmd)
    
  }
  
}


## bioscience
for (i in 1:length(met$A.CEL_file))
{
  print(paste(met$A.CEL_file[i],date()))
  
  ##open the cufflink file
  if(file.exists(paste0("/scratch2/brain_express/tophatJuly13/Sample_",met$A.CEL_file[i],"/trim_galore/align_summary.txt")))
  {  
    cmd <- paste0("grep 'Mapped' /scratch2/brain_express/tophatJuly13/Sample_",met$A.CEL_file[i],"/trim_galore/align_summary.txt |  cut -f17 -d' ' | awk '{ sum+=$1} END {print sum}'  >> ~/test.txt")
    system(cmd )
    
  }
  
}



















