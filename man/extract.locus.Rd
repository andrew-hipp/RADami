\name{extract.locus}
\alias{extract.locus}
\title{Extract sequence data for one or a few loci}
\description{
  This function is pretty rudimentary and I've used it but little. As much as anything, it's a reminder for how to get sequence data out of a \code{pyRAD.loci} object.
}
\usage{
extract.locus(pyDat, locName, dat.format = c("text", "fasta", "matrix"))
}
\arguments{
  \item{pyDat}{an object of class \code{pyRAD.loci}}
  \item{locName}{character name of the locus or loci you want}
  \item{dat.format}{currently takes text, fasta, or matrix, but the fasta format would still need to be written to a file}
}
\details{}
\value{a named character vector if \code{dat.format %in% c('text', 'fasta'), otherwise a matrix}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrew Hipp
}
