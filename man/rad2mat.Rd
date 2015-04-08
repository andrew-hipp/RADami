\name{rad2mat}
\alias{rad2mat}
\title{Generates a matrix of sequences from a \code{pyRAD.loci} object}
\description{This function exports a matrix of sequences from a \code{pyRAD.loci} object,
			 with aligned sequences in the cells, individuals as rows, and loci as columns.
			 As a matrix, this is very easy to work with in R for subsetting and creating 
			 new datasets.}
\usage{
rad2mat(pyDat, fill.N = TRUE)
}
\arguments{
  \item{pyDat}{an object of class \code{pyRAD.loci}}
  \item{fill.N}{indicates whether empty cells should be filled with Ns (\code{TRUE}) or not}
}
\details{In general, filling with Ns is a handy thing to do, even though it presumably eats up
		 a bunch of memory. As written, RADami is not memory-friendly.}
\value{A matrix of sequences from a \code{pyRAD.loci} object,
			 with aligned sequences in the cells, individuals as rows, and loci as columns}
\author{Andrew Hipp}
\seealso{
	\code{\link{rad2phy}},
	\code{\link{summary.pyRAD.loci}}
}
\examples{
data(oak.rads)
out <- rad2mat(oak.rads)
}