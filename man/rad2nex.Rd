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
  \item{fillBlanks}{integer; number of characters to fill blanks with}
  \item{fillChar}{character to fill blanks with}
  \item{loci}{loci to include in the Nexus file, by name}
  \item{outfile}{file name to write}
  \item{verbose}{boolean; report number of characters to add to file}
  \item{logfile}{name of logfile, which holds the names of loci in the phylip file written}
}
\value{No value. This function is called for the side-effect of writing a Nexus file}
\author{Andrew Hipp}
\seealso{
	\code{\link{rad2mat}},
  \code{\link{rad2phy}},
	\code{\link{write.DNAStringSet}}
}
\keyword{IO}
