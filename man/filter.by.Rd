\name{filter.by}
\alias{filter.by}
\title{Identify names of loci for which at the requested taxa are present,
       either all or some number of them}
\description{Returns the names of loci for which the requested taxa are present at some threshold}
\usage{
filter.by(dat, taxa, req,
          threshInds = 2, threshTaxa = 2,
          use.tidyName = FALSE, ...)
}
\arguments{
  \item{dat}{an object of class \code{pyRAD.loci} or \code{summary.pyRAD.loci}}
  \item{taxa}{a list of character vectors of taxa or individuals}
  \item{req}{a boolean vector corresponding to taxa, indicating whether the taxon is required (TRUE) or optional (FALSE)}
  \item{threshInds}{a vector of integers corresponding to taxa, indicating the threshold number of individuals}
  \item{threshTaxa}{a single integer, indicating how many taxa are required at their \code{threshInds} level to export a locus}
  \item{use.tidyName}{whether to use the \code{tidyName} function for matching taxon names;
                    useful if there is a mix of uppercase, lowercase, special characters etc.
                    that are not really essential for distinguishing taxa}
  \item{\dots}{extra arguments to \code{tidyName}}
}
\details{Not currently called by any other functions in this package, but useful for quickly
         identifying loci for which data are present in a specified number of individuals}
\value{A character vector of locus names}
\author{Andrew Hipp}
\seealso{
  \code{\link{subset.pyRAD.loci}}
}
\keyword{manip}
