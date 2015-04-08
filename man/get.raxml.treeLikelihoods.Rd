\name{get.raxml.treeLikelihoods}
\alias{get.raxml.treeLikelihoods}
\title{Extract likelihoods from a RAxML info file}
\description{Reads the info file from a RAxML site-likelihood analysis with multiple input trees.
			 Probably not often needed on its own, but used in \code{match.lnL.to.trees}.}
\usage{
get.raxml.treeLikelihoods(x, logfile = NA)
}
\arguments{
  \item{x}{file name of a RAxML .info file from site-likelihood analysis}
  \item{logfile}{name of a log file, useful for recording any files that were not successfully
				 read in}
}
\value{A named vector of class \code{double} with tree likelihoods, where the names are character 
       equivalents of the tree numbers; or, if the file had no trees in it, the character vector
	   "FAIL"}
\author{Andrew Hipp}
\seealso{
	\code{\link{match.lnL.to.trees}},
	\code{\link{get.raxml.siteLikelihoods}}
}
\keyword{IO}