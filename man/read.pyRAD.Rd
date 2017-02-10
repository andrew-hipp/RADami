\name{read.pyRAD}
\alias{read.pyRAD}
\title{Read in data from pyRAD}
\description{This function begins essentially all analyses in this package. It reads in the
			 pyRAD loci file and stores it in a form to facilitate a variety of downstream
			 analyses coded in this package.}
\usage{
read.pyRAD(filename, reportInterval = 20000,
	         breakLinesSeparate = FALSE, 
           doSummary = TRUE, ...)
}
\arguments{
  \item{filename}{name of the pyRAD loci file}
  \item{reportInterval}{how often to report progress while reading and processing files}
  \item{breakLinesSeparate}{are there separate break lines (blank) and summary lines (not blank)
							between loci in the pyRAD file?}
  \item{doSummary}{should \code{summary.pyRAD.loci} be run on the object and appended to the end?}
  \item{\dots}{extra arguments to be passed along to \code{readLines}}
}
\details{As written, this function has a lot of summary embedded in it. The default \code{doSummary = TRUE}
		 should be left for the other higher-level functions to work properly}
\value{
	A list with the following items:
		\item{tips}{a character vector of the tips for each sequence}
		\item{seqs}{a character vector of all sequences}
		\item{breaks}{an integer vector indexing the lines that are breaks between loci}
		\item{break.vectors}{a named character vector of the break lines between vectors, with names
							 corresponding to locus names}
		\item{cons}{a named character vector of the lines summarizing each locus, with names
					corresponding to locus names; if \code{breakLinesSeparate = F}, this is identical to
					\code{break.vectors}}
		\item{locus.index}{an index corresponding to \code{seqs} and \code{tips}, indicating
						   what locus each belongs to}
		\item{file.read}{name of the file read in}
		\item{timestamp}{timestamp for exiting the function}
		\item{radSummary}{an object of class \code{summary.pyRAD.loci}, if called, which it is by default}
}
\author{Andrew Hipp}
\seealso{
	\code{\link{summary.pyRAD.loci}}
}
\keyword{IO}
