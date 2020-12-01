\name{rad2nex}
\alias{rad2nex}
\title{Export a \code{rad.mat} object to a Nexus file}
\description{
  This function takes output from \code{rad2mat} and writes it to a
  Nexus file, ready for analysis by PAUP*, BEAST, or other software
  }
\usage{
rad2nex(pyDat, inds = row.names(pyDat), indNames = NA,
        fillBlanks = NA, fillChar = '-',
        loci = dimnames(pyDat)[[2]],
        outfile = 'pyMat.out.nex',
        verbose = FALSE,
        logfile = 'rad2nex.log', ...)
}
\arguments{
  \item{pyDat}{an object of class \code{rad.mat}}
  \item{inds}{individuals to include in the export}
  \item{indNames}{
    a vector of names to rename individuals if desired;
    assumed to be in same order as inds
    }
  \item{fillBlanks}{...}
  \item{fillChar}{...}
  \item{loci}{...}
  \item{outfile}{...}
  \item{verbose}{...}
  \item{logfile}{name of logfile, which holds the names of loci in the phylip file written}
}
\value{No value. This function is called for the side-effect of writing a phylip file}
\author{Andrew Hipp}
\seealso{
	\code{\link{rad2mat}},
	\code{\link{write.DNAStringSet}}
}
\keyword{IO}
