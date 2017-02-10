\name{subset.pyRAD.loci}
\alias{subset.pyRAD.loci}
\title{Subsets a \code{pyRAD.loci} object}
\description{This function returns a dataset (at this point, only in
             \code{DNAStringSet} format) from a \code{pyRAD.loci} object,
			 subsetting by loci and taxa provided as arguments}
\usage{
\method{subset}{pyRAD.loci}(x, loci = colnames(x$radSummary$inds.mat),
              taxa = row.names(x$radSummary$inds.mat),
		format = "DNAStringSet", reportInterval = 500, mins = 1,
		nucVarType = c("verystrict", "strict", "relaxed"),
		use.tidyName = FALSE,
		snpsOnly = FALSE, cores = 1, ...)
}
\arguments{
  \item{x}{an object of class \code{pyRAD.loci}}
  \item{loci}{a character vector of loci to include; defaults to all loci}
  \item{taxa}{a character vector of taxa to include; defaults to all taxa}
  \item{format}{format for DNA matrix; at this point, only \code{DNAStringSet}
				is supported}
  \item{reportInterval}{interval (number of loci) for reporting progress}
  \item{mins}{minimum number of individuals required to return locus data}
  \item{nucVarType}{nucleotide variability to consider in screening for polymorphism:
					'relaxed' for all nucleotides, including ambiguities
					other than 'N', and counting all variable sites; 'strict'
					for only A,C,G,T, still counting all variable sites; and
					'verystrict' for only A,C,G,T, countin only parsimony-
					informative sites.}
  \item{use.tidyName}{choose whether names are passed through tidyName}
  \item{snpsOnly}{choose whether you only want to get SNPs back, or whole locus sequences}
  \item{cores}{how many processing cores to use; if \code{cores = 1}, \code{lapply} is used
                throughout instead of \code{mclapply}}
  \item{\dots}{additional arguments passed to tidyName}
}
\details{Probably it would make sense to wrap some of the functioning of
		 \code{gen.RAD.loci.datasets} into this function, but for now they
		 are separated.}
\value{
	A list with values:
		\item{DNA}{a list of datasets, one per locus, defined by \code{format}
				   argument}
		\item{variable}{a boolean vector, indicating which loci are variable according
						to the criterion provided in arguments}
		\item{ntaxa}{an integer vector, indicating how many taxa are in each locus}
}
\author{Andrew Hipp}
\seealso{
	\code{\link{gen.RAD.loci.datasets}},
	\code{\link{filter.by}}
	}
