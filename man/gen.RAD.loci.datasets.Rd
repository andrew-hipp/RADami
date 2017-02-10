\name{gen.RAD.loci.datasets}
\alias{gen.RAD.loci.datasets}
\title{Generate datasets for partitioned RAD analysis}
\description{This function takes a set of trees (typically from \code{genTrees}), RAD loci,
             and some parameters about how taxa to include in the datasets, how many taxa
			 at minimum need to be included in the datasets, and whether only variable loci
			 will be included.}
\usage{
gen.RAD.loci.datasets(rads, trees = "none", loci = "all",
                      taxa = "all", minTaxa = 4, 
                      onlyVariable = TRUE, fileBase = "DEFAULT",
                      splitInto = 1,
                      cores = 2,
                      raxSinglePath = "raxmlHPC-AVX",
                      raxMultiPath = "raxmlHPC-PTHREADS-AVX",
                      header = "#!/bin/sh")
}
\arguments{
  \item{rads}{an object of class \code{pyRAD.loci}}
  \item{trees}{a set of trees, or 'none' if no trees are being used (e.g., if you just
               want to export loci one-by-one for analysis)}
  \item{loci}{a vector of locus names to include in analysis, or 'all.' Typically, use 'all'
              here, because part of what this function does is subset loci for you based
			  on criteria you provide}
  \item{taxa}{a vector of taxa to include in analysis, or 'all' for all taxa}
  \item{minTaxa}{an integer, specifying the minimum number of taxa required to export a locus}
  \item{onlyVariable}{either \code{TRUE} or \code{FALSE} depending on whether loci should be
                      variable or not to be exported}
  \item{fileBase}{base name for the directories and files created by this function;
                  defaults to a handy dated format}
  \item{splitInto}{the number of folders to split files into for analysis, in case you
                   want to run several instances of RAxML to get through analysis more
				   quickly}
  \item{cores}{the number of processing cores to use for the RAxML full-matrix analysis}
  \item{raxSinglePath}{the path to your RAxML executable, single-core installation}
  \item{raxMultiPath}{the path to your RAxML executable, multi-core installation}
  \item{header}{header for the top of your shell scripts / bash files; not necessary
                for a Windows batch file (\code{.bat}), but defaults to the header
				needed for a Linux shell script}
}
\details{This function subsets the treeset you give it for all the loci imported, then
         matched loci and treesets and a batch file for calculating the likelihood of all
         trees for each locus. It is called for the side-effect of writing all these files.

         All files are exported to a series of folders equal in number to \code{splitInto}
		 + 1. As currently written, the folder with the suffix \code{.0} holds all batch files, an index of trees
		 that matches the trees you hand to the function (as \code{trees}) with the trees
		 it returns in each subset, and the \code{.info} output files that are used by
		 \code{get.raxml.treeLikelihoods} as well as the site-likelihood files. The
		 alignments and trees, as well as the reduced alignments created by RAxML, remain
		 in the split folders (1 through \code{splitInto}).

		 This function is the workhorse of the partitioned RAD analysis and is time-consuming to run.
         It could be executed in parallel and should be in future versions. What is most time-consuming
		 is that for every treeset exported, the function makes sure that all trees are unique
		 (see \code{unique.multiPhylo} from \code{ape}).}
\value{no value. Called for the side effect of the files it exports.}
\references{Hipp A.L., Eaton D.A.R., Cavender-Bares J., Fitzek E., Nipper R. & Manos P.S. (Accepted pending revision). A framework phylogeny of the American oak clade based on sequenced RAD data. PLoS ONE.}
\author{Andrew Hipp}
\seealso{
  \code{\link{genTrees}},
  \code{\link{get.raxml.treeLikelihoods}},
  \code{\link{match.lnL.to.trees}},
  \code{\link{rank.partitionedRAD}}
  }
\keyword{datagen}
\keyword{manip}
