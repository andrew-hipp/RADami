\name{do.a.rad.comparison.tree}
\alias{do.a.rad.comparison.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
do.a.rad.comparison.tree(tr, rads, cex.scalar = 2, analysis = "desc.og", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{
%%     ~~Describe \code{tr} here~~
}
  \item{rads}{
%%     ~~Describe \code{rads} here~~
}
  \item{cex.scalar}{
%%     ~~Describe \code{cex.scalar} here~~
}
  \item{analysis}{
%%     ~~Describe \code{analysis} here~~
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
function (tr, rads, cex.scalar = 2, analysis = "desc.og", ...) 
{
    nodes.comparison <- do.all.nodes(tr, rads, analysis)
    plot(tr)
    nodelabels(pie = nodes.comparison[, "prop"], node = as.numeric(row.names(nodes.comparison)), 
        cex = cex.scalar * nodes.comparison[, "abs"]/max(nodes.comparison[, 
            "abs"]), piecol = c("black", "white"), ...)
    return(nodes.comparison)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
