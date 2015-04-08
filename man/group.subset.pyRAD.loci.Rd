\name{group.subset.pyRAD.loci}
\alias{group.subset.pyRAD.loci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
group.subset.pyRAD.loci(dat, groups, mins = 10, loci = names(dat$DNA), use.tidyName = TRUE, cores = 8, ...)
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
function (dat, groups, mins = 10, loci = names(dat$DNA), use.tidyName = TRUE, 
    cores = 8, ...) 
{
    if (use.tidyName) {
        if (cores > 1) 
            dat$DNA <- mclapply(dat$DNA, function(x) {
                names(x) <- tidyName(names(x), ...)
                x
            }, mc.cores = cores)
        else dat$DNA <- lapply(dat$DNA, function(x) {
            names(x) <- tidyName(names(x), ...)
            x
        })
        groups <- lapply(groups, tidyName, ...)
    }
    if (cores > 1) 
        out <- t(mcmapply(function(x) c(sapply(groups, function(y) sum(names(x) \%in\% 
            y)), total = length(x)), dat$DNA, mc.cores = cores))
    else out <- t(sapply(dat$DNA, function(x) c(sapply(groups, 
        function(y) sum(names(x) \%in\% y)), total = length(x))))
    if (!is.na(mins[1])) 
        leave.in <- apply(out, 1, function(x) all(x >= mins))
    else leave.in <- rep(TRUE, dim(out)[1])
    out <- out[leave.in, ]
    if (!is.na(loci[1])) 
        out <- out[intersect(row.names(out), loci), ]
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
