\name{rtreePhylo}
\alias{rtreePhylo}
\title{Generate random trees using an input tree to get taxon names and tips}
\description{A wrapper around \code{rtree}}
\usage{
rtreePhylo(x, N = 1, ...)
}
\arguments{
  \item{x}{an object of class \code{phylo}}
  \item{N}{number of trees to return}
  \item{\dots}{extra arguments to pass along to \code{rtree}}
}
\value{an object of class \code{multiPhylo}}
\author{Andrew Hipp}
\seealso{
	\code{\link{genTrees}}
}
\examples{
	data(oak.tree)
	a <- rtreePhylo(oak.tree, 4)
	layout(matrix(1:4,2,2))
	lapply(a, plot)
	}