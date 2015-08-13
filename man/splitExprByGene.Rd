\name{splitExprByGene}
\alias{splitExprByGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splitExprByGene(i, ensemblRef, ensemblGenes, PUTM, exprPUTM, SNIG, exprSNIG, pathResExpr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{ensemblRef}{
%%     ~~Describe \code{ensemblRef} here~~
}
  \item{ensemblGenes}{
%%     ~~Describe \code{ensemblGenes} here~~
}
  \item{PUTM}{
%%     ~~Describe \code{PUTM} here~~
}
  \item{exprPUTM}{
%%     ~~Describe \code{exprPUTM} here~~
}
  \item{SNIG}{
%%     ~~Describe \code{SNIG} here~~
}
  \item{exprSNIG}{
%%     ~~Describe \code{exprSNIG} here~~
}
  \item{pathResExpr}{
%%     ~~Describe \code{pathResExpr} here~~
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
function (i, ensemblRef, ensemblGenes, PUTM, exprPUTM, SNIG, 
    exprSNIG, pathResExpr) 
{
    sys.source("/home/adai/scripts/common_functions.R", attach(NULL, 
        name = "myenv"))
    ensemblRef.gene <- ensemblGenes[i, ]
    rm(i)
    geneID <- ensemblRef.gene$Ensembl.Gene.ID
    gene.info <- ensemblRef[which(ensemblRef$Ensembl.Gene.ID == 
        geneID), ]
    print(geneID)
    dir.create(pathResExpr, showWarnings = FALSE)
    fn.rda <- paste0(pathResExpr, "/byGene_snps1Mb/", geneID, 
        ".rda")
    dir.create(dirname(fn.rda), showWarnings = FALSE)
    tissues <- c("PUTM", "SNIG")
    expr <- vector(mode = "list", length(tissues))
    names(expr) <- tissues
    IDs <- PUTM[which(PUTM$A.CEL_file \%in\% as.character(rownames(exprPUTM))), 
        1]
    IDs <- gsub("/", "_", IDs)
    rownames(exprPUTM) <- IDs
    rm(indID, IDs)
    IDs <- SNIG[which(SNIG$A.CEL_file \%in\% as.character(rownames(exprSNIG))), 
        1]
    IDs <- gsub("/", "_", IDs)
    rownames(exprSNIG) <- IDs
    rm(IDs)
    expr[["SNIG"]] <- as.matrix(exprSNIG[, is.element(colnames(exprSNIG), 
        ensemblRef.gene$Ensembl.Gene.ID)])
    expr[["PUTM"]] <- as.matrix(exprPUTM[, is.element(colnames(exprPUTM), 
        ensemblRef.gene$Ensembl.Gene.ID)])
    if (length(expr[["SNIG"]]) == 0 && length(expr[["PUTM"]]) == 
        0) {
        cat(geneID, "\n", file = paste0(pathResExpr, "/genesNoExpressed"), 
            append = T)
        stop("Gene not expressed")
    }
    for (tissue in tissues) {
        if (length(expr[[tissue]]) == 0) {
            cat(geneID, "\n", file = paste0(pathResExpr, "/genesNoExpressed", 
                tissue), append = T)
            tissues <- tissues[-which(tissues == tissue)]
        }
        if (ncol(expr[["SNIG"]]) == 1) {
            colnames(expr[["SNIG"]]) <- intersect(colnames(exprSNIG), 
                ensemblRef.gene$Ensembl.Gene.ID)
        }
        if (ncol(expr[["PUTM"]]) == 1) {
            colnames(expr[["PUTM"]]) <- intersect(colnames(exprPUTM), 
                ensemblRef.gene$Ensembl.Gene.ID)
        }
    }
    save(expr, gene.info, file = fn.rda, compress = "bzip2", 
        ascii = T)
    rm(fn.rda)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
