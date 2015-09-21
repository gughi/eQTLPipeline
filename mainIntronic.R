###  Manuel Sebastian Guelfi
### 21-5-15
### This script qunatify the abundance of the intronic region in each gene.

### the approach we will use is to load the expression quantification of the exonic regions only 
### produced by HTSeq-counts  and then we will load the expression of the genes accounting for
### exons and introns (quantified with intersect bed) and the we will subtract the first dataset.

## load the reads for exonic regions
load("data/expr/rawCounts/genic/exprSQ.rda")
## dim(SQcounts)
## [1]   170 64078

## load the reads for intronic and exonic regions
SQEIcounts <- read.csv("data/expr/rawCounts/genic/exprExonIntr.csv",row.names=1)
## dim(SQEIcounts)
## [1]   170 64078

SQEIcounts <- SQEIcounts[as.character(rownames(SQcounts)),]

## we check that the order of the column and rownames have the same order
stopifnot(identical(colnames(SQcounts),colnames(SQEIcounts)))
stopifnot(identical(rownames(SQcounts),rownames(SQEIcounts)))

## we calculate the intronic reads
intronicReads <- SQEIcounts - SQcounts
## dim(intronicReads)
## [1]   170 64078

summedInt <- apply(intronicReads,1,sum)
summedIntExo <- apply(SQEIcounts,1,sum)
summedExo <- apply(SQcounts,1,sum)

ratiosInt <- summedInt/summedIntExo
ratiosExo <- summedExo/summedIntExo

jpeg("/home/guelfi/tmp/test.jpeg", type="cairo")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(ratiosInt,col='skyblue',border=F,main= "Proportion of genic reads" ,xlab="" ,xaxt="n")
barplot(ratiosExo,add=T,col='red',border=F,xlab= "Samples",xaxt="n")
legend("topright", inset=c(-0.3,0), legend=c("Intronic","Exonic"), 
       fill=c("skyblue","red"), col=c("skyblue","red"), title="Reads type")
dev.off()
