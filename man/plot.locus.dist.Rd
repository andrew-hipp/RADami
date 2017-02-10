\name{plot.locus.dist}
\alias{plot.locus.dist}
\title{Plot a phylogeny with a graphical representation of locus sharing}
\description{Plots a phylogeny lined up with the locus-sharing similarity matrix in graphical
			 style}
\usage{
\method{plot}{locus.dist}(x, tr, trW = 3, plotW = 5, labelsW = 3, plotGap = 0.25,
     scalar = 1.5, barH = 1, barW = 10, bar.cex = 0.5,
		 point.pch = c(21, 21),
     cols = c('black', 'red'), ...)
}
\arguments{
  \item{x}{output from \code{locus.dist}}
  \item{tr}{an object of class \code{phylo} with tips named the same as the rows of \code{locD}}
  \item{trW}{width of the panels holding the tree}
  \item{plotW}{width of the panel holding the matrix plot}
  \item{labelsW}{width of the panel holding the taxon labels (right side of plot)}
  \item{plotGap}{gap between the side panels and the central panel}
  \item{scalar}{the scalar for the dots in the locus distance plot}
  \item{barH}{height of the panel holding the bars}
	\item{barW}{width of individual bars; goes to \code{lwd}}
	\item{bar.cex}{font size for numbers above the bars}
  \item{point.pch}{plot character (\code{pch}) for the off-diagonal and diagonal elements of the matrix plot}
  \item{cols}{color for the off-diagonal and diagonal elements of the matrix plot}
  \item{\dots}{additional arguments passed along to \code{text}}
}
\details{
	Use default options in \code{locus.dist} to get a sensible plot.
	Also, while different plotting characters can be defined for the off-diagonal and diagonal
	elements of the matrix plot, these may be hard to compare with each other in size; probably
	the diagonal elements are better distinguished by color or shading.

	Currently, this function is implemented using \code{layout}, which makes it rather clunky
	and inflexible.
	}
\value{
	A list comprising:
		\item{Xs}{the x-coordinates of dots in the matrix plot}
		\item{heights}{heights of the bars at the top of the plot, which show the overlap of
					   each individual averaged over all other individuals on the plot}

}
\author{Andrew Hipp}
\seealso{
	\code{\link{plot.rankedPartitionedRAD}},
	\code{\link{locus.dist}}
}
\examples{
  data(oak.rads)
  data(oak.tree)
  plot.locus.dist(locus.dist(oak.rads, proportional = TRUE), oak.tree)
  # and now a nicer-looking plot:
  plot.locus.dist(locus.dist(oak.rads, proportional = TRUE, upper = TRUE),
                  ladderize(oak.tree), trW = 1, labelsW = 1, scalar = 3)
}
\keyword{hplot}
