\name{group.pyRAD.loci}
\alias{group.pyRAD.loci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
group.pyRAD.loci(dat, groups, mins = 10, loci = dimnames(dat$radSummary$inds.mat)[[2]], use.tidyName = TRUE, cores = 8, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{mins}{
%%     ~~Describe \code{mins} here~~
}
  \item{loci}{
%%     ~~Describe \code{loci} here~~
}
  \item{use.tidyName}{
%%     ~~Describe \code{use.tidyName} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (dat, groups, mins = 10, loci = dimnames(dat$radSummary$inds.mat)[[2]], 
    use.tidyName = TRUE, cores = 8, ...) 
{
    if (use.tidyName) {
        row.names(dat$radSummary$inds.mat) <- tidyName(row.names(dat$radSummary$inds.mat), 
            ...)
        groups <- lapply(groups, tidyName, ...)
    }
    inds.mats <- mclapply(groups, function(x) temp <- dat$radSummary$inds.mat[x[x \%in\% 
        row.names(dat$radSummary$inds.mat)], ], mc.cores = cores)
    out <- mcmapply(function(x) apply(x, 2, sum), inds.mats, 
        mc.cores = cores)
    if (!is.na(mins[1])) 
        leave.in <- apply(out, 1, function(x) all(x >= mins))
    else leave.in <- rep(TRUE, dim(out)[1])
    out <- cbind(out[leave.in, ], total = colSums(dat$radSummary$inds.mat)[leave.in])
    if (!is.na(loci[1])) 
        out <- out[intersect(row.names(out), loci), ]
    attr(out, "groups") = groups
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
