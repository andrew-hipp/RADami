\name{match.lnL.to.trees}
\alias{match.lnL.to.trees}
\title{Match the likelihoods returned on a subset of trees for numerous loci to the trees originally input}
\description{
			 This function reads in tree likelihoods for a series of loci exported from
			 \code{gen.RAD.loci.datasets} and matches those up to the trees originally created.
			 Because \code{gen.RAD.loci.datasets} only exports trees that include the taxa in
			 each locus exported, these likelihoods are matched up to all trees in the original
			 dataset that are compatible with each pruned tree.
			}
\usage{
match.lnL.to.trees(directory = "getwd()",
                   lnLprefix = "RAxML_info.",
		   lnLsuffix = ".lnL",
                   treeIndexFile = "tree.index.lines.txt",
		   locus.names = NULL, ...)
}
\arguments{
  \item{directory}{path for the directory in which the 'info' and 'lnL' files were saved by RAxML}
  \item{lnLprefix}{prefix for the RAxML info files; generally should not be changed, unless RAxML
                   changes this down the line}
  \item{lnLsuffix}{suffix appended to the files; this is currently should not be changed, as it is
				   hard-coded in \code{gen.RAD.loci.datasets}}
  \item{treeIndexFile}{names of the tree index file; again, this file name is hard-coded in
					   \code{gen.RAD.loci.datasets}}
  \item{locus.names}{a vector of locus names, if you want to subset by loci at the outset of
					 analysis; typically, subsetting is just done by the function}
  \item{\dots}{additional parameters, not currently passed along to anybody}
}
\details{The idea behind this function is that we want to let loci vote only on trees
		 they are qualified to vote on. In an earlier implementation, each locus was forced to
		 vote on every tree, unpruned, and the results were simply noisy. Here, we take a set of
		 trees and, for each locus out a treeset pruned to just that taxon set for each taxon
		 set. We then calculate the site likelihoods based on those treesets.

		 Every locus is, in fact, evaluating every tree, but because many trees become
		 topologically identical after pruning, we have to set an equal score to trees that
		 are topologically identical... and which trees are topologically identical differs
		 among loci.}
\value{A matrix of likelihoods, with loci as rows and the original trees as columns}
\references{Hipp A.L., Eaton D.A.R., Cavender-Bares J., Fitzek E., Nipper R. & Manos P.S. (Accepted pending revision). A framework phylogeny of the American oak clade based on sequenced RAD data. PLoS ONE.}
\author{Andrew Hipp}
\seealso{
	\code{\link{gen.RAD.loci.datasets}},
	\code{\link{plot.rankedPartitionedRAD}},
	\code{\link{rank.partitionedRAD}}
}
\keyword{manip}
