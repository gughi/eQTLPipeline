\name{doSwamp}
\alias{doSwamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doSwamp(RPKM.cqn, covs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RPKM.cqn}{
%%     ~~Describe \code{RPKM.cqn} here~~
}
  \item{covs}{
%%     ~~Describe \code{covs} here~~
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
function (RPKM.cqn, covs) 
{
    library(swamp)
    RPKM.cqn.tmp <- t(RPKM.cqn)
    expr.data.o <- RPKM.cqn.tmp[order(rownames(RPKM.cqn.tmp)), 
        ]
    traits.o <- covs[as.character(rownames(expr.data.o)), ]
    rownames(traits.o) <- rownames(expr.data.o)
    print(head(covs))
    stopifnot(identical(rownames(traits.o), rownames(expr.data.o)))
    res1 <- prince(t(expr.data.o), traits.o, top = 15, permute = TRUE)
    prince.plot(prince = res1)
    res2 <- prince.var.plot(t(expr.data.o), show.top = 50, npermute = 10)
    hca.plot(as.matrix(t(expr.data.o)), as.data.frame(traits.o))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
