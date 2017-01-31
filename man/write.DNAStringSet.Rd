\name{write.DNAStringSet}
\alias{write.DNAStringSet}
\title{Writes a sequence matrix to phylip format}
\description{
This function writes a Biostrings \code{DNAStringSet} object to a file. Currently, phylip and fasta formats are supported
}
\usage{
write.DNAStringSet(x, format= c('phylip', 'fasta'), padding = 30,
                  filename = "DNAStringSetOut.phy",
                  fastaPrefix = ">")
}
\arguments{
  \item{x}{an object of class \code{DNAStringSet}}
  \item{format}{at this point, only fasta and phylip formats are supported}
  \item{padding}{total number of characters, including sequence title, to place before each sequence}
  \item{filename}{full name of file and path}
  \item{fastaPrefix}{prefix character to add to taxon label for fasta files}
}
\value{At this point, returns a 0 to indicate successful completion}
\author{Andrew Hipp}
\seealso{
  \code{\link{rad2phy}}
}
