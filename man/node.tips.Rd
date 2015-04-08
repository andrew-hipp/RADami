\name{node.tips}
\alias{node.tips}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
node.tips(tr, node, getOut = c("desc.og", "desc"), add.og = ">AC_h")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{
%%     ~~Describe \code{tr} here~~
}
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{getOut}{
%%     ~~Describe \code{getOut} here~~
}
  \item{add.og}{
%%     ~~Describe \code{add.og} here~~
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
function (tr, node, getOut = c("desc.og", "desc"), add.og = ">AC_h") 
{
    to.do <- tr$edge[which(tr$edge[, 1] == node), 2]
    out <- lapply(Descendants(tr, to.do, "tips"), function(x) tr$tip.label[x])
    if (getOut[1] == "desc.og") 
        out[[3]] <- c(tr$tip.label[!tr$tip.label \%in\% unique(unlist(out))], 
            add.og)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
