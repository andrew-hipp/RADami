\name{compare.loci}
\alias{compare.loci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compare.loci(taxa, radMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxa}{
%%     ~~Describe \code{taxa} here~~
}
  \item{radMat}{
%%     ~~Describe \code{radMat} here~~
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
function (taxa, radMat) 
{
    if (length(taxa[[1]]) > 1) 
        v1 <- apply(radMat[taxa[[1]], ], 2, any)
    else v1 <- radMat[taxa[[1]], ]
    if (length(taxa[[2]]) > 1) 
        v2 <- apply(radMat[taxa[[2]], ], 2, any)
    else v2 <- radMat[taxa[[2]], ]
    if (length(taxa) == 3) {
        if (length(taxa[[3]]) > 1) 
            v3 <- apply(radMat[taxa[[3]], ], 2, any)
        else v3 <- radMat[taxa[[2]], ]
        v.mat <- do.call(rbind, list(v1, v2, v3))
    }
    else v.mat <- rbind(v1, v2)
    v.out <- sum(apply(v.mat, 2, all))
    v.denom <- sum(apply(v.mat, 2, any))
    out <- c(abs = v.out, prop = v.out/v.denom)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
