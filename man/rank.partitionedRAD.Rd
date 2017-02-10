\name{rank.partitionedRAD}
\alias{rank.partitionedRAD}
\title{Bin trees into 'suppported' or 'disfavored' by locus.}
\description{This function takes output from \code{match.lnL.to.trees} and scores trees as
			 favored or disfavored for all loci that satisfy the criteria specified}
\usage{
rank.partitionedRAD(radMat, minTrees = 10,
	                  min.overall.diff.lnL = 5,
		threshold.lnL = 2,
                    discardDoubleCounts = TRUE)
}
\arguments{
  \item{radMat}{matrix output from \code{match.lnL.to.trees}}
  \item{minTrees}{integer, minimum number of unique trees required per locus}
  \item{min.overall.diff.lnL}{minimum log-likelihood difference required between the most poorly-
							  supported tree and the best supported tree for each locus}
  \item{threshold.lnL}{the log-likelihood window for binning a tree into favored or disfavored}
  \item{discardDoubleCounts}{determines whether that place any tree within both the favored and
							 disfavored bins should be counted (\code{TRUE}) or not
							 (\code{FALSE})}
}
\details{This function works by first filtering loci by the minimum number of trees and the
		 overall range in log-likelihood, then binning trees into supported if they are within
		 \code{threshold.lnL} of the best-supported tree for each locus, and disfavored if they
		 are within \code{threshold.lnL} of the least supported tree for each locus. Overlapping
		 trees can occur when the range in log-likelihood is set to less than twice the threshold
		 log-likelihood or quite near it.}
\value{
	A list composed of 7 matrices and one numeric vector. The matrices all have loci as the rows and trees as the columns:
		\item{bestMat}{boolean; which trees are supported by which loci}
		\item{worstMat}{boolean; which trees are disfavored by which loci}
		\item{doubleCountMat}{boolean; which trees are double-counted by which loci}
		\item{radMat}{the final likelihood matrix corresponding with \code{bestMat} and \code{worstMat}}
		\item{radMat.preLnLDiff}{the likelihood matrix before filtering by lnL range}
		\item{radMat.preMinTrees}{the likelihood matrix before filtering by minimum number of trees}
		\item{radMat.preDoubleCounts }{the likelihood matrix before filtering by double counts}
		\item{params}{the analysis parameters, in this order:
					  minTrees, min.overall.diff.lnL, threshold.lnL, discardDoubleCounts}
}
\references{Hipp A.L., Eaton D.A.R., Cavender-Bares J., Fitzek E., Nipper R. & Manos P.S. (Accepted pending revision). A framework phylogeny of the American oak clade based on sequenced RAD data. PLoS ONE.}
\author{Andrew Hipp}
\seealso{
	\code{\link{match.lnL.to.trees}},
	\code{\link{plot.rankedPartitionedRAD}}
}
