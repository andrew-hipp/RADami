\name{genotypes.pyRAD.loci}
\alias{genotypes.pyRAD.loci}
\title{Generates genotypes from a \code{subset.pyRAD.loci} object}
\description{Makes a dataframe of SNP calls from a pyRAD.loci object for export to hierfstat}
\usage{
genotypes.pyRAD.loci(dat, groups, loci = 'all', taxa = 'all',
	                             useSnps = c('first', 'all'), concat = c(FALSE, TRUE),
					use.tidyName = TRUE, na.rm = c('none', 'columns', 'rows'), maxAlleles = 2,
					tidyVals = c('-', '.','>', '_', ' '), sortByGroups = TRUE,
					variable.only = FALSE, make.dummy.column = TRUE, alleleDigits = 2, toInteger = TRUE, missingData = "0000",
					cores = 1)
}
\arguments{
  \item{dat}{an object of class \code{subset.pyRAD.loci}}
  \item{groups}{a list of vectors, each vector designating a set of individuals to count for each group}
  \item{loci}{a character vector of loci to include; defaults to all loci}
  \item{taxa}{a character vector of taxa to include; defaults to all taxa}
  \item{useSnps}{xxx}
  \item{concat}{xxx}
  \item{use.tidyName}{choose whether or not to use \code{tidyName} when matching
                       names}
  \item{na.rm}{xxx}
  \item{maxAlleles}{xxx}
  \item{tidyVals}{xxx}
  \item{sortByGroups}{xxx}
  \item{variableOnly}{xxx}
  \item{make.dummy.column}{xxx}
  \item{alleleDigits}{xxx}
  \item{toInteger}{xxx}
  \item{missingData}{xxx}
  \item{cores}{how many processing cores to use; if \code{cores = 1}, \code{lapply} is used
                throughout instead of \code{mclapply}}
}
\details{current not tightly integrated with hierfstat, but designed to be used by it}
\value{
	data.frame of genotypes
}
\author{Andrew Hipp}
\seealso{
	\code{\link{subset.pyRAD.loci}}
	}
