\name{overlap.report}
\alias{overlap.report}
\title{Generates a report on locus overlap between pairs of samples named in a structured fashion}
\description{Provided a \code{pyRAD.loci} object with a subset of paired sample names that differ 
			 only by a suffix (defined here as \code{repPattern} and \code{origPattern}), 
			 this function will return a table summarizing the overlap pattern}
\usage{
overlap.report(dat, repPattern = "_re", origPattern = "_h")
}
\arguments{
  \item{dat}{an object of class \code{pyRAD.loci} or \code{summary.pyRAD.loci}, or
             just the \code{indsMat} object from a \code{summary.pyRAD.loci} object}
  \item{repPattern}{suffix for the first half of each pair}
  \item{origPattern}{suffix for the second half of each pair}
}
\details{This may not be of use to too many people, but we used it in the article below. In
		 our case, the samples were paired technical replicates, but you could imagine sister
		 species, populations or anything else being the replicates of interest.}
\value{a matrix with pairs of individuals as rows and overlap statistics as columns}
\references{Hipp A.L., Eaton D.A.R., Cavender-Bares J., Fitzek E., Nipper R. & Manos P.S. (Accepted pending revision). A framework phylogeny of the American oak clade based on sequenced RAD data. PLoS ONE.}
\author{Andrew Hipp}
\seealso{
	\code{\link{help}}, ~~~
}
\examples{
	data(oak.rads)
	overlap.report(oak.rads)
}