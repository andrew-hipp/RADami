\name{genTrees}
\alias{genTrees}
\title{
Generate permuted trees using rNNI or rtree, and check for uniqueness.
}
\description{
This function generates a specified number of trees by permutation or random
generation, then checks to make sure that they are unique before returning the
original tree plus the generated trees
}
\usage{
genTrees(x, N = 200, filebase = "trial", method = c("nni", "random"),
         maxmoves = 2, perms = "DEFAULT", software = c("raxml"), ...)
}
\arguments{
  \item{x}{an object of class \code{"phylo"}}
  \item{N}{number of new trees to return}
  \item{filebase}{base name of the file to which trees will be written}
  \item{method}{method of generating the trees, either nearest neighbor interchange (nni) or random}
  \item{maxmoves}{maximum number nni moves}
  \item{perms}{number of permutations per nni stratum; defaults to the number of single-step nni trees,
               then the remainder of N divided by the number of nni strata requested (negatives not allowed)}
  \item{software}{not currently used; in an earlier version, would write out analysis file for PAUP or RAxML}
  \item{\dots}{further arguments to be passed to rtree}
}
\details{
The main need for this function is to ensure that the trees you get back are unique... mainly just a wrapper to rNNI and rtree.
\code{perms} defaults to:

\code{c(length(nni(x)), max(1, as.integer(N - (length(nni(x))))))}
}
\value{an object of class \code{"multiPhylo"}}
\author{Andrew Hipp}
\seealso{
  \code{\link{gen.RAD.loci.datasets}},
  \code{\link{compare.all.trees}}
}
\examples{
  data(oak.tree)
  out <- genTrees(oak.tree, N = 10, maxmoves = 2, perms = c(5, 5))
  plot(out[[1]]) # this is the tree you input
  plot(out[[2]]) # this is the first generated tree...
  plot(out[[11]]) # ... and this is the last generated tree.
}
