correPlot <-
function (mat1,mat2,xlab)
{
  library(gplots)
  linp<-matrix(ncol=ncol(mat1),nrow=ncol(mat2))  
  rownames(linp)<-colnames(mat2)
  colnames(linp)<-colnames(mat1)
  rsquared<-matrix(ncol=ncol(mat1),nrow=ncol(mat2))
  rownames(rsquared)<-colnames(mat2)
  colnames(rsquared)<-colnames(mat1)
  for (i in 1:ncol(mat2)){
    for (j in 1:ncol(mat1)){
      fit<-lm(mat1[,j]~mat2[,i])
      s<-summary(fit)
      linp[i,j]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
      rsquared[i,j]<-s$r.squared[1]
    }}
  
  
  smallest=-20
  linp10<-log10(linp)
  linp10<-replace(linp10,linp10<=smallest,smallest)
  print(linp10)
  
  rsquaredTMP <- rsquared
  rsquaredTMP[which(linp10>-5)] = NA
  rsquaredTMP <-round(rsquaredTMP,digits=2)
  heatmap.2(linp10,Colv=F,Rowv=F,dendrogram="none",trace="none",symbreaks=F,symkey=F,breaks=seq(-20,0,length.out=100),key=T,col=heat.colors(99),
            cexRow=1,cexCol=1,colsep=NULL,rowsep=NULL,sepcolor=sepcolor,sepwidth=sepwidth,
            main="",labCol=paste(1:ncol(linp10),sep=""),margins=c(5,7),labRow=,xlab=xlab,
            cellnote=rsquaredTMP,notecol="black",notecex=1) 
#   heatmap.2(linp10,Colv=F,Rowv=F,dendrogram="none",trace="none",symbreaks=F,symkey=F,breaks=seq(-20,0,length.out=100),key=T,col=heat.colors(99),
#             cexRow=1,cexCol=1,colsep=NULL,rowsep=NULL,sepcolor=sepcolor,sepwidth=sepwidth,
#             main="",labCol=paste(1:ncol(linp10),sep=""),margins=c(5,7),labRow=,xlab=xlab,
#             cellnote=matrix(ncol=ncol(linp),nrow=nrow(linp)),notecol="black",notecex=1) 
#   
  
  
}
