\name{plot.rankedPartitionedRAD}
\alias{plot.rankedPartitionedRAD}
\title{Plot a \code{rankedPartitionedRAD} object}
\description{This function creates a plot with 1-3 panels
			 describing the number of loci supporting, disfavoring,
			 or ambiguous for each tree, plotted against log-likelihood}
\usage{
\method{plot}{rankedPartitionedRAD}(x, 
     tree.lnL.file = NULL,
     fileprefix = NULL, 
     lnL.break = NULL,
     regression = NULL,
     ci = NULL,
     widthScalar = 0.85, 
     panels = c("bestMat", "worstMat", "doubleCountMat"), 
     squareSize = switch(as.character(length(panels)), '2' = 5, '3' = 3), 
     primeTreeColor = "red", 
     primeTreeCharacter = 19, 
     highlight.points = NULL, 
     highlight.colors = NULL,
     filebase = 'DEFAULT',
     ...)
}
\arguments{
  \item{x}{an object of class \code{rankedPartitionedRAD}}
  \item{tree.lnL.file}{file name of the info file written by RAxML; optional}
  \item{fileprefix}{prefix of the files to save; as written, you have to provide some prefix
                    to identify the files if you want a pdf, but the remainder of the file name 
					defaults to a combination of parameters and date}
  \item{lnL.break}{a numeric vector, which should be the length of \code{panels};
				   if not null, the lowest tree likelihood to show on each panel}
  \item{regression}{a boolean vector, which should be the length of \code{panels};
				   if not null, whether to draw a linear fit line on each panel}
  \item{ci}{a numeric vector, the length of \code{panels}; 
            what CI to show on each regression, or \code{0} for no CI}
  \item{widthScalar}{a multiplier on the width of the entire plot to get the plotting panels square
					 or not as you prefer}
  \item{panels}{a character vector describing the order of panels; 
				can be any of "bestMat", "worstMat", "doubleCountMat" in any order}
  \item{squareSize}{the size of the panels in inches; defaults to 5" for a 2-panel plot, 
					3" for a 3-panel plot, only because that works well in my trials}
  \item{primeTreeColor}{the color of the dot representing the first tree in the matrix, which 
						should be your original tree (unpermuted) if you've been following the 
						workflow from \code{genTrees} through this final and exciting plot}
  \item{primeTreeCharacter}{the character of the dot representing the first tree in the matrix}
  \item{highlight.points}{number of any trees to highlight}
  \item{highlight.colors}{colors vector or single color for highlighted points}
  \item{filebase}{the base of the file name for exported plots}
  \item{\dots}{additional arguments passed to \code{plot}}
}
\details{The x-axis will either be the tree likelihoods estimated by RAxML in a separate analysis (e.g., using the entire dataset),
		 or the summed locus likelihoods from the partitioned RAD analysis if no filename is provided.
		 This function as written prints to pdf if you enter a \code{fileprefix}, to the screen if you don't.
		 Default filebase pastes \code{rad.partitioned} with the date and the parameter values used. 
		 If default filebase is not desired, enter something other than "DEFAULT."}
\value{This function is called primarily for its side-effect of printing plots to a file.
		It returns, invisibly, a list of the fit objects for each panel, only relevant 
		if you set any elements of \code{regression} to \code{TRUE}.}
\references{Hipp A.L., Eaton D.A.R., Cavender-Bares J., Fitzek E., Nipper R. & Manos P.S. (Accepted pending revision). A framework phylogeny of the American oak clade based on sequenced RAD data. PLoS ONE.}
\author{Andrew Hipp}
\seealso{
	\code{\link{rank.partitionedRAD}},
	\code{\link{plot.locus.dist}}
}
\keyword{hplot}