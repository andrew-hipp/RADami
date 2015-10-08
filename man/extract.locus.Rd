\name{extract.locus}
\alias{extract.locus}
\title{Extract sequence data for one or a few loci}
\description{
  This function is pretty rudimentary and I've used it but little. As much as anything, it's a reminder for how to get sequence data out of a \code{pyRAD.loci} object.
}
\usage{
extract.locus(pyDat, locName, dat.format = c("text", "fasta", "matrix"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pyDat}{an object of class \code{pyRAD.loci}}
  \item{locName}{character name of the locus or loci you want}
  \item{dat.format}{currently takes text, fasta, or matrix, but the fasta format would still need to be written to a file}
}
\details{}
\value{a named character vector if \code{dat.format %in% c('text', 'fasta'), otherwise a matrix}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrew Hipp
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
