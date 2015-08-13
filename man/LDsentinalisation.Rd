\name{LDsentinalisation}
\alias{LDsentinalisation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LDsentinalisation(i, ensemblGenes, pathFinalSentinalised, pathUnsentinalised, FDRthr, exprLocation, snpLocation, genotypeFile, tmpFolder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{ensemblGenes}{
%%     ~~Describe \code{ensemblGenes} here~~
}
  \item{pathFinalSentinalised}{
%%     ~~Describe \code{pathFinalSentinalised} here~~
}
  \item{pathUnsentinalised}{
%%     ~~Describe \code{pathUnsentinalised} here~~
}
  \item{FDRthr}{
%%     ~~Describe \code{FDRthr} here~~
}
  \item{exprLocation}{
%%     ~~Describe \code{exprLocation} here~~
}
  \item{snpLocation}{
%%     ~~Describe \code{snpLocation} here~~
}
  \item{genotypeFile}{
%%     ~~Describe \code{genotypeFile} here~~
}
  \item{tmpFolder}{
%%     ~~Describe \code{tmpFolder} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (i, ensemblGenes, pathFinalSentinalised, pathUnsentinalised, 
    FDRthr, exprLocation, snpLocation, genotypeFile, tmpFolder) 
{
    sys.source("/home/adai/scripts/common_functions.R", attach(NULL, 
        name = "myenv"))
    geneID <- ensemblGenes[i, 1]
    rm(i)
    print(geneID)
    dir.create(file.path(pathFinalSentinalised), showWarnings = FALSE)
    dir.create(file.path(pathFinalSentinalised, "SNIG"), showWarnings = FALSE)
    dir.create(file.path(pathFinalSentinalised, "PUTM"), showWarnings = FALSE)
    fn.rda <- paste0(exprLocation, geneID, ".rda")
    load(fn.rda)
    fn.rda <- paste0(snpLocation, geneID, ".rda")
    load(fn.rda)
    rm(fn.rda)
    tissues <- c("SNIG", "PUTM")
    library(MatrixEQTL)
    for (tissue in tissues) {
        print(tissue)
        if (!file.exists(paste0(pathUnsentinalised, tissue, "/", 
            geneID))) {
            next
        }
        eQTLs <- read.delim(paste0(pathUnsentinalised, tissue, 
            "/", geneID), sep = " ")
        j = 1
        my.cov0 <- NA
        while (TRUE) {
            if (j > 1) {
                if (!file.exists(paste0(tmpFolder, tissue, "/", 
                  geneID))) {
                  break
                }
                eQTLs <- read.delim(paste0(tmpFolder, tissue, 
                  "/", geneID), sep = " ")
            }
            else {
                system(paste("echo", "snps", "gene", "t-stat", 
                  "pvalue", "FDR", "beta", "myFDR", "degree", 
                  paste0(">> ", pathFinalSentinalised, tissue, 
                    "/", geneID)))
            }
            if (min(eQTLs$myFDR) > FDRthr) {
                print(min(eQTLs$myFDR))
                break
            }
            if (file.exists(paste0(tmpFolder, tissue, "/", geneID))) {
                system(paste0("rm ", tmpFolder, tissue, "/", 
                  geneID))
            }
            SNP <- eQTLs$snps[which(eQTLs$myFDR == min(eQTLs$myFDR))[1]]
            leadEQTL <- eQTLs[which(eQTLs$myFDR == min(eQTLs$myFDR))[1], 
                ]
            print(eQTLs$snps[which(eQTLs$myFDR == min(eQTLs$myFDR))])
            tmpRes <- data.frame(cbind(eQTLs[which(eQTLs$myFDR == 
                min(eQTLs$myFDR)), ], j))
            for (n in 1:nrow(tmpRes)) {
                system(paste("echo", tmpRes$snps[n], tmpRes$gene[n], 
                  tmpRes$statistic[n], tmpRes$pvalue[n], tmpRes$FDR[n], 
                  tmpRes$beta[n], tmpRes$myFDR[n], tmpRes$j[n], 
                  paste0(">> ", pathFinalSentinalised, tissue, 
                    "/", geneID)))
            }
            my.expr0 <- as.matrix(expr[[tissue]])
            my.expr0 <- as.matrix(my.expr0[, which(colMeans(is.na(my.expr0)) != 
                1)])
            if (ncol(my.expr0) < 1) {
                next
            }
            my.markers0 <- as.matrix(markers[, is.element(colnames(markers), 
                rownames(my.expr0))])
            if (j < 2) {
                my.covTMP <- read.table.rows(paste0(genotypeFile), 
                  keepRows = rownames(my.expr0), sep = " ", header = F)
                my.cov0 <- as.matrix(my.covTMP[colnames(my.markers0), 
                  2:4])
                my.cov0 <- t(my.cov0)
            }
            covDose <- markers[which(SNP == rownames(markers)), 
                colnames(my.cov0)]
            stopifnot(identical(colnames(covDose), colnames(my.cov0)))
            my.cov0 <- rbind(my.cov0, covDose)
            my.expr0 <- as.matrix(my.expr0[colnames(my.markers0), 
                ])
            if (ncol(my.expr0) == 1) {
                colnames(my.expr0) <- colnames(expr[[tissue]])
            }
            my.expr0 <- t(my.expr0)
            stopifnot(identical(colnames(my.expr0), colnames(my.markers0)))
            stopifnot(identical(colnames(my.cov0), colnames(my.markers0)))
            my.expr <- SlicedData$new()
            my.expr$CreateFromMatrix(my.expr0)
            my.markers <- SlicedData$new()
            my.markers$CreateFromMatrix(my.markers0)
            my.cov <- SlicedData$new()
            my.cov$CreateFromMatrix(as.matrix(my.cov0))
            rm(my.expr0, my.markers0)
            outputFile = paste0(tmpFolder, tissue, "/", geneID)
            store <- Matrix_eQTL_main(my.markers, my.expr, my.cov, 
                output_file_name = NULL, pvOutputThreshold = 0.1, 
                useModel = modelLINEAR, errorCovariance = numeric(0), 
                verbose = T)
            pval <- store$all$eqtls$pvalue
            myFDR <- sapply(pval, function(pval) p.adjust(pval, 
                method = "fdr", n = nrow(markers.info)))
            rm(pval)
            my.eQTLstmp <- cbind(store$all$eqtls, myFDR)
            eQTLs <- my.eQTLstmp[which(my.eQTLstmp$myFDR <= FDRthr), 
                ]
            rm(my.eQTLstmp)
            if (nrow(eQTLs) > 0) {
                write.table(eQTLs, outputFile, row.names = F)
            }
            else {
                break
            }
            rm(eQTLs)
            j = j + 1
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
