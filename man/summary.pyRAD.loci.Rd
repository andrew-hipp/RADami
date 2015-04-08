\name{summary.pyRAD.loci}
\alias{summary.pyRAD.loci}
\title{Summary method for \code{pyRAD.loci} objects}
\description{This function provides the summary stats used by many functions in \code{RADami}}
\usage{
\method{summary}{pyRAD.loci}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{pyRAD.loci}}
  \item{\dots}{not currently used, but maintained for consistency with default \code{summary} method}
}
\details{This is currently called by default by \code{read.pyRAD}, creating a list object called \code{radSummary}
		 in value for that object.}
\value{
	A list with of the following objects:
		\item{num.loci}{an integer, indicating how many loci there are in the object}
		\item{tips.per.locus}{a list of character vectors, enumerating the tips in each locus}
		\item{break.vectors}{a named character vector of the break lines between vectors, with names
							 corresponding to locus names -- completely redundant with the same object in the 
							 \code{pyRAD.loci} object}
		\item{seqs.per.locus}{a list of named character vectors, enumerating the sequences in each locus, named
							  by tips}
		\item{num.inds.per.locus}{an integer vector, indicating how many individuals are in each locus}
		\item{variable.loci}{a boolean vector, indicating which loci are variable (only using the pyRAD
							 break vectors)}
		\item{inds.mat}{a handy boolean matrix with individuals as the rows and loci as the columns; 
						used in several places in \code{RADami}}
		\item{locus.lengths}{an integer vector, indicating the length of each locus}
}
\author{Andrew Hipp}
\note{Probably these functions could be teased apart to be more useful and save space, but for now
		 they are left as is. At this point, there is actually a fair bit of redundancy.}
\seealso{
	\code{\link{read.pyRAD}},
	\code{\link{print.pyRAD.loci}}
}