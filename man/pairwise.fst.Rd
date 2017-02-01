\name{pairwise.fst.loci}
\alias{pairwise.fst.loci}
\title{Returns pairwise FST list for groups of individuals}
\description{for use with RADami}
\usage{
pairwise.fst.loci(dat, group.list, to.do,
                  minInds = 3, do.hist = TRUE, cores = 1)
}
\arguments{
  \item{dat}{an object of class \code{subset.pyRAD.loci}}
  \item{group.list}{named list of character vectors, indicating groups of individuals}
  \item{to.do}{which groups to generate pairwise data for}
  \item{do.hist}{boolean; whether a histogram of pairwise FST should be written to the working directory}
  \item{cores}{how many processing cores to use}
}
\value{
A list with the following elements:
  \item{fst}{pairwise fst list}
  \item{groups}{groups used in the analysis}
  \item{genotypes}{genotype data}
  \item{sp.pairs}{species pairs analyzed}
  \item{snpLocs}{location of SNPs}
  \item{timestamp}{when analysis was conducted}
}
\author{Andrew Hipp}
\seealso{
	\code{\link{subset.pyRAD.loci}}
}
