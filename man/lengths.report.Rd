\name{lengths.report}
\alias{lengths.report}
\title{Return the lengths of all loci in a pyRAD.loci object}
\description{Return the lengths of all loci in a pyRAD.loci object}
\usage{
lengths.report(dat, numtodo = 10,
               reportInterval = 2000,
               high.mem = TRUE)
}
\arguments{
  \item{dat}{an object of class \code{pyRAD.loci}}
  \item{numtodo}{an integer specifying how many loci to return the length of}
  \item{reportInterval}{how many loci between reports of time required to complete when using
                        \code{high.mem} = \code{FALSE}}
  \item{high.mem}{forces the code to loop through instead of using \code{sapply}}
}
\details{Called by \code{summary.pyRAD.loci}, and thus not ordinarily needed by users.
		 If looping is used to get lengths, the function runs very slowly, and
		 \code{reportInterval} can be useful. In early tests, I was having memory allocation
		 issues, so \code{high.mem} was useful. I haven't needed it since then, but
		 maybe you will.
		 }
\value{A named vector of class \code{integer}, with names corresponding to loci whose lengths are returned}
\author{Andrew Hipp}
\seealso{
	\code{\link{read.pyRAD}},
	\code{\link{summary.pyRAD.loci}}
}
\examples{
  data(oak.rads)
  lengths.report(oak.rads, -1) # returns all loci
  lengths.report(oak.rads, 20) # returns only the first 20 loci
  identical(oak.rads$radSummary$locus.lengths, lengths.report(oak.rads, -1))
}
