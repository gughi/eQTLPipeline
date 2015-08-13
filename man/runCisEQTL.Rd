\name{runCisEQTL}
\alias{runCisEQTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runCisEQTL(i, ensemblGenes, exprLocation, snpLocation, outputFolder, genotypeFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{ensemblGenes}{
%%     ~~Describe \code{ensemblGenes} here~~
}
  \item{exprLocation}{
%%     ~~Describe \code{exprLocation} here~~
}
  \item{snpLocation}{
%%     ~~Describe \code{snpLocation} here~~
}
  \item{outputFolder}{
%%     ~~Describe \code{outputFolder} here~~
}
  \item{genotypeFile}{
%%     ~~Describe \code{genotypeFile} here~~
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
function (i, ensemblGenes, exprLocation, snpLocation, outputFolder, 
    genotypeFile) 
{
    sys.source("/home/adai/scripts/common_functions.R", attach(NULL, 
        name = "myenv"))
    geneID <- ensemblGenes[i, 1]
    rm(i)
    geneID <- sub(".rda", "", geneID)
    print(geneID)
    fn.rda <- paste0(exprLocation, geneID, ".rda")
    load(fn.rda)
    fn.rda <- paste0(snpLocation, geneID, ".rda")
    load(fn.rda)
    rm(fn.rda, gene.info, t.map)
    dir.create(file.path(outputFolder), showWarnings = FALSE)
    dir.create(file.path(outputFolder, "resMatrixEQTL/"), showWarnings = FALSE)
    dir.create(file.path(paste0(outputFolder, "resMatrixEQTL/"), 
        "SNIG"), showWarnings = FALSE)
    dir.create(file.path(paste0(outputFolder, "resMatrixEQTL/"), 
        "PUTM"), showWarnings = FALSE)
    tissues <- c("PUTM", "SNIG")
    if (markers == "No polymorphisms found within +/- 1Mb of TSS or TES") {
        stop("No polymorphisms found within +/- 1Mb of TSS or TES")
    }
    library(MatrixEQTL)
    for (tissue in tissues) {
        print(tissue)
        my.expr0 <- as.matrix(expr[[tissue]])
        my.expr0 <- as.matrix(my.expr0[, which(colMeans(is.na(my.expr0)) != 
            1)])
        if (ncol(my.expr0) < 1) {
            next
        }
        my.markers0 <- as.matrix(markers[, is.element(colnames(markers), 
            rownames(my.expr0))])
        my.covTMP <- read.table.rows(paste0(genotypeFile), keepRows = rownames(my.expr0), 
            sep = " ", header = F)
        my.cov0 <- as.matrix(my.covTMP[colnames(my.markers0), 
            2:4])
        my.cov0 <- t(my.cov0)
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
        my.cov$CreateFromMatrix(my.cov0)
        rm(my.expr0, my.markers0, my.cov0)
        outputFile = paste0(outputFolder, "resMatrixEQTL/", tissue, 
            "/", geneID)
        store <- Matrix_eQTL_main(my.markers, my.expr, my.cov, 
            output_file_name = NULL, pvOutputThreshold = 0.1, 
            useModel = modelLINEAR, errorCovariance = numeric(0), 
            verbose = T)
        pval <- store$all$eqtls$pvalue
        myFDR <- sapply(pval, function(pval) p.adjust(pval, method = "fdr", 
            n = nrow(markers.info)))
        rm(pval)
        my.eQTLstmp <- cbind(store$all$eqtls, myFDR)
        my.eQTLs <- my.eQTLstmp[which(my.eQTLstmp$myFDR <= 0.1), 
            ]
        rm(my.eQTLstmp)
        print(nrow(my.eQTLs))
        if (nrow(my.eQTLs) > 0) {
            write.table(my.eQTLs, outputFile, row.names = F)
        }
    }
    rm(my.expr, my.markers, my.cov, markers.info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
