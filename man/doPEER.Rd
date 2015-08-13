\name{doPEER}
\alias{doPEER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doPEER(RPKM.cqn, nFactors, covs, outputFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RPKM.cqn}{
%%     ~~Describe \code{RPKM.cqn} here~~
}
  \item{nFactors}{
%%     ~~Describe \code{nFactors} here~~
}
  \item{covs}{
%%     ~~Describe \code{covs} here~~
}
  \item{outputFile}{
%%     ~~Describe \code{outputFile} here~~
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
function (RPKM.cqn, nFactors, covs, outputFile) 
{
    library(peer, lib = "/home/seb/Rlibrary/")
    covs <- covs[as.vector(rownames(RPKM.cqn)), ]
    stopifnot(identical(rownames(covs), rownames(RPKM.cqn)))
    model = PEER()
    PEER_setNk(model, nFactors)
    PEER_setCovariates(model, data.matrix(covs))
    PEER_setPhenoMean(model, as.matrix(RPKM.cqn))
    print(paste("dimension of the expression matrix", print(dim(PEER_getPhenoMean(model)))))
    print(paste("dimension of the covariate matrix", print(dim(PEER_getCovariates(model)))))
    PEER_update(model)
    factor <- PEER_getX(model)
    head(factor)
    rownames(factor) <- rownames(RPKM.cqn)
    write.csv(factor, file = outputFile)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
