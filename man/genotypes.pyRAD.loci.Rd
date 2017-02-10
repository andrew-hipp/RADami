\name{genotypes.pyRAD.loci}
\alias{genotypes.pyRAD.loci}
\title{Generates genotypes from a \code{subset.pyRAD.loci} object}
\description{Makes a dataframe of SNP calls from a pyRAD.loci object for export to hierfstat}
\usage{
genotypes.pyRAD.loci(dat, groups, loci = 'all', taxa = 'all',
	        useSnps = c('first', 'all'), concat = c(FALSE, TRUE),
					use.tidyName = FALSE,
					na.rm = c('none', 'columns', 'rows'),
					maxAlleles = 2,
					tidyVals = c('-', '.','>', '_', ' '),
					sortByGroups = TRUE,
					variable.only = FALSE, make.dummy.column = TRUE,
					alleleDigits = 1,
					toInteger = TRUE, missingData = "00",
					cores = 1)
}
\arguments{
  \item{dat}{an object of class \code{subset.pyRAD.loci}}
  \item{groups}{a list of vectors, each vector designating a set of individuals to count for each group}
  \item{loci}{a character vector of loci to include; defaults to all loci}
  \item{taxa}{a character vector of taxa to include; defaults to all taxa}
  \item{useSnps}{whether to use the first SNP in each locus or all}
  \item{concat}{boolean: whether to concatenate SNPs or return by locus}
  \item{use.tidyName}{choose whether or not to use \code{tidyName} when matching names}
  \item{na.rm}{whether to remove missin data by columns, by rows, or not at all}
  \item{maxAlleles}{maximum number of alleles per locus}
  \item{tidyVals}{characters to be stripped out by tidyName, if using}
  \item{sortByGroups}{boolean: whether to sort snps by groups}
  \item{variable.only}{boolean: whether to keep only variable sites}
  \item{make.dummy.column}{boolean: whether to add a dummy column when there is only one SNP
	                         (currently added as a work-around for \code{wc})}
  \item{alleleDigits}{how many digits per allele}
  \item{toInteger}{boolean: whether to return SNPs as integers}
  \item{missingData}{string to use for missing data}
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
