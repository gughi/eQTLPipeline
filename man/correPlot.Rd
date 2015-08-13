\name{correPlot}
\alias{correPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correPlot(mat1, mat2, xlab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat1}{
%%     ~~Describe \code{mat1} here~~
}
  \item{mat2}{
%%     ~~Describe \code{mat2} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
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
function (mat1, mat2, xlab) 
{
    library(gplots)
    linp <- matrix(ncol = ncol(mat1), nrow = ncol(mat2))
    rownames(linp) <- colnames(mat2)
    colnames(linp) <- colnames(mat1)
    rsquared <- matrix(ncol = ncol(mat1), nrow = ncol(mat2))
    rownames(rsquared) <- colnames(mat2)
    colnames(rsquared) <- colnames(mat1)
    for (i in 1:ncol(mat2)) {
        for (j in 1:ncol(mat1)) {
            fit <- lm(mat1[, j] ~ mat2[, i])
            s <- summary(fit)
            linp[i, j] <- pf(s$fstatistic[1], s$fstatistic[2], 
                s$fstatistic[3], lower.tail = FALSE)
            rsquared[i, j] <- s$r.squared[1]
        }
    }
    smallest = -20
    linp10 <- log10(linp)
    linp10 <- replace(linp10, linp10 <= smallest, smallest)
    print(linp10)
    rsquaredTMP <- rsquared
    rsquaredTMP[which(linp10 > -5)] = NA
    rsquaredTMP <- round(rsquaredTMP, digits = 2)
    heatmap.2(linp10, Colv = F, Rowv = F, dendrogram = "none", 
        trace = "none", symbreaks = F, symkey = F, breaks = seq(-20, 
            0, length.out = 100), key = T, col = heat.colors(99), 
        cexRow = 1, cexCol = 1, colsep = NULL, rowsep = NULL, 
        sepcolor = sepcolor, sepwidth = sepwidth, main = "", 
        labCol = paste(1:ncol(linp10), sep = ""), margins = c(5, 
            7), labRow = , xlab = xlab, cellnote = rsquaredTMP, 
        notecol = "black", notecex = 1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
