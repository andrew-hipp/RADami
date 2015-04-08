\name{write.DNAStringSet}
\alias{write.DNAStringSet}
\title{Writes a sequence matrix to phylip format}
\description{
This function writes a Biostrings \code{DNAStringSet} object to a file. Currently, only phylip format is supported
}
\usage{
write.DNAStringSet(x, format = "phylip", padding = 30, filename = "DNAStringSetOut.phy")
}
\arguments{
  \item{x}{an object of class \code{DNAStringSet}}
  \item{format}{at this point, only phylip format is supported}
  \item{padding}{total number of characters, including sequence title, to place before each sequence}
  \item{filename}{full name of file and path}
}
\value{At this point, returns a 0 to indicate successful completion}
\author{Andrew Hipp}
\seealso{
  \code{\link{rad2phy}}
}
