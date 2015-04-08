\name{print.pyRAD.loci}
\alias{print.pyRAD.loci}
\title{Print method for \code{pyRAD.loci} objects}
\description{Prints very basic summary data for a \code{pyRAD.loci} object 
			 rather than streaming the entire thing, which can be huge and not 
			 so useful to read on the screen}
\usage{
	\method{print}{pyRAD.loci}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{pyRAD.loci}}
  \item{...}{unused}
}
\value{No value. Called for the side-effect of printing attractively to the screen}
\author{Andrew Hipp}
\seealso{
  \code{\link{read.pyRAD}}
}