\name{extract.locus}
\alias{extract.locus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract.locus(pyDat, locName, dat.format = c("text", "fasta", "matrix"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pyDat}{
%%     ~~Describe \code{pyDat} here~~
}
  \item{locName}{
%%     ~~Describe \code{locName} here~~
}
  \item{dat.format}{
%%     ~~Describe \code{dat.format} here~~
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
function (pyDat, locName, dat.format = c("text", "fasta", "matrix")) 
{
    seqs <- pyDat$seqs[pyDat$locus.index == locName]
    tips <- pyDat$tips[pyDat$locus.index == locName]
    if (dat.format[1] == "text") 
        names(seqs) <- pyDat$tips[pyDat$locus.index == locName]
    if (dat.format[1] == "fasta") 
        seqs <- as.character(matrix(c(tips, seqs), 2, byrow = T))
    if (dat.format[1] == "matrix") 
        seqs <- as.matrix(seqs)
    seqs
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
