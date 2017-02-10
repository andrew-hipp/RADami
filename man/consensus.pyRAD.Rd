\name{consensus.pyRAD}
\alias{consensus.pyRAD}
\title{Generate a consensus sequence for each of the loci in a pyRAD.loci object}
\description{Uses \code{consensusString} [S4 method for signature 'DNAStringSet']
             from the \code{Biostrings} package to create
             a consensus sequence for each locus in a \code{pyRAD.loci} object. All options
			 in \code{consensusString} can be passed in. As a side-effect, writes
			 a fasta-style file to \code{writeFile}.}
\usage{
consensus.pyRAD(pyIn, from = NA, to = NA, fastaNames = T,
    writeFile = 'rads.con.txt', cores = 1, ...)
}
\arguments{
  \item{pyIn}{an object of class \code{pyRAD.loci}}
  \item{from}{first locus to summarize (by number, sequentially from the beginning of the file),
              or \code{NA} if all loci are to be used}
  \item{to}{last locus to summarize (by number, sequentially from the beginning of the file),
            or \code{NA} if all loci are to be used}
  \item{fastaNames}{use \code{TRUE} to append a \code{>} to the front of each sequence title, making a fasta-
                    style file}
  \item{writeFile}{filename to which to write the sequences}
  \item{cores}{number of cores to use; requires the \code{parallel} package}
  \item{\dots}{extra arguments to \code{consensusString}}
}
\value{A character vector of sequences, with locus names as the vector names.}
\author{Andrew Hipp}
\examples{
  data(oak.rads)
  # provide filename for export in line below
  temp = consensus.pyRAD(oak.rads, from=1, to = 10, writeFile = NA)
  temp
}
\keyword{manip}
