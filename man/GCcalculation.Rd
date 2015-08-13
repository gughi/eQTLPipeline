\name{GCcalculation}
\alias{GCcalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GCcalculation(region, genRef)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{region}{
%%     ~~Describe \code{region} here~~
}
  \item{genRef}{
%%     ~~Describe \code{genRef} here~~
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
function (region, genRef) 
{
    print("creting the BED file")
    tmpBED <- tempfile("GCcont", fileext = ".BED")
    BED <- paste(gsub("chr", "", region[, 1]), region[, 2], region[, 
        3], rownames(region), sep = "\t")
    write.table(data.frame(BED), file = tmpBED, row.names = F, 
        col.names = F, quote = F)
    rm(BED)
    tmpGCcon <- tempfile("GCcont")
    cmd <- paste0("bedtools nuc -fi ", genRef, " -bed ", tmpBED, 
        " > ", tmpGCcon)
    print("executing the bedtools")
    system(cmd)
    rm(cmd)
    print("collecting the GC content")
    GCcontent <- read.delim(pipe(paste("cut -f4,6", tmpGCcon)))
    colnames(GCcontent) <- c("ID", "GCcontent")
    rownames(GCcontent) <- GCcontent$ID
    GCcontent$ID <- NULL
    GCcontent
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
