\name{rad2phy}
\alias{rad2phy}
\title{Export a \code{rad.mat} object to a phylip file}
\description{This function takes output from \code{rad2mat} and writes it to a phylip file, ready for analysis by RAxML or other software}
\usage{
rad2phy(pyDat, inds = row.names(pyDat), loci = dimnames(pyDat)[[2]],
        outfile = "pyMat.out.phy", padding = 50, verbose = FALSE,
        logfile ="rad2phy.log")
}
\arguments{
  \item{pyDat}{an object of class \code{rad.mat}}
  \item{inds}{individuals to include in the export}
  \item{loci}{loci to include in the export}
  \item{outfile}{file name for export}
  \item{padding}{total number of characters, including sequence title, to place before each sequence}
  \item{verbose}{whether to report on loci being written}
  \item{logfile}{name of logfile, which holds the names of loci in the phylip file written}
}
\value{No value. This function is called for the side-effect of writing a phylip file}
\author{Andrew Hipp}
\seealso{
	\code{\link{rad2mat}},
	\code{\link{write.DNAStringSet}}
}
\keyword{IO}
