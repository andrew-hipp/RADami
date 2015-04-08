\name{get.raxml.siteLikelihoods}
\alias{get.raxml.siteLikelihoods}
\title{
Read likelihoods from a RAxML_perSiteLLs file
}
\description{
This function is used to import likelihoods from a RAxML site likelihoods export. 
}
\usage{
get.raxml.siteLikelihoods(x)
}
\arguments{
  \item{x}{Filename of the RAxML_perSiteLLs file}
}
\details{
Not used in any downstream analyses currently in RADbuddy.
}
\value{A matrix with trees as rows, site lnL as columns}
\author{Andrew Hipp}
\seealso{
  \code{\link{get.raxml.treeLikelihoods}},
  \code{\link{match.lnL.to.trees}} 
}
\keyword{manip}
