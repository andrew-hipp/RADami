\name{locus.dist}
\alias{locus.dist}
\title{Generate a pairwise similarity matrix for individuals based on locus-sharing patterns}
\description{This function returns a pairwise similarity matrix based on how many loci or
			 the proportion of loci shared between individuals.}
\usage{
locus.dist(pyIn, proportional = TRUE, upper = TRUE, diagonal = TRUE)
}
\arguments{
  \item{pyIn}{an object of class \code{pyRAD.loci}}
  \item{proportional}{set to \code{TRUE} for a proportional similarity (scaled from 0 to 1),
					  \code{FALSE} for absolute number of loci shared}
  \item{upper}{set to \code{TRUE} for a symmetrical matrix, \code{FALSE} for a lower-only matrix
			   }
  \item{diagonal}{set to \code{TRUE} to get the number (or percentage) of total loci reported
				  for each individual}
}
\details{Under the \code{proportional = TRUE} option, this function normalizes based on all loci
		 in the dataset. Thus, two individuals that share 100% of their loci, but both of which
		 have sequences for only 90% of the loci reported across the entire dataset, will report
		 0.9 instead of 1.0.}
\value{A square matrix with individuals as the rows and columns, and proportions or whole numbers
	   in the cells}
\author{Andrew Hipp}
\seealso{
	\code{\link{plot.locus.dist}}
}
\examples{
	data(oak.rads)
	temp <- locus.dist(oak.rads, proportional = TRUE)
	diag(temp) # this is how well everyone sequenced as a proportion
	temp <- locus.dist(oak.rads, proportional = FALSE)
	temp[c(1,28),]# how well did a successful sequencing do in comparison to a poor sequencing?
	temp[6:7,] # ... and 2010 sequencing (upper row) vs. 2012 sequencing of the same library?
}
