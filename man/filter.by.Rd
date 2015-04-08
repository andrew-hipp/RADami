\name{filter.by}
\alias{filter.by}
\title{Identify names of loci for which at the requested taxa are present, 
       either all or some number of them}
\description{Returns the names of loci for which the requested taxa are present at some threshold}
\usage{
filter.by(dat, taxa, threshold = "all")
}
\arguments{
  \item{dat}{an object of class \code{pyRAD.loci} or \code{summary.pyRAD.loci}}
  \item{taxa}{a character vector of taxa or individuals to identify in \code{dat}}
  \item{threshold}{an integer setting the minimum number of taxa needed for each locus, 
                   or "all" for all taxa}
}
\details{Not currently called by any other functions in this package, but useful for quickly
         identifying loci for which data are present in a specified number of individuals}
\value{A character vector of locus names}
\author{Andrew Hipp}
\seealso{ 
  \code{\link{subset.pyRAD.loci}}
}
\keyword{manip}
