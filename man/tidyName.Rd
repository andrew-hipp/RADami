\name{tidyName}
\alias{tidyName}
\title{
Clean up a name
}
\description{
Deletes in a character vector all the characters provided in \code{fixes}.
}
\usage{
tidyName(x, fixes = c("_", ".", " ", "-", ">"), case = c("lower", "upper", "nochange"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the vector in which to make replacements
}
  \item{fixes}{
all characters to be replaced
}
  \item{case}{
change case to upper, lower, or make no changes
}
}
\details{
This function currently uses no regular expressions or wildcards and is made deliberately very thorough, tedious even, as it cycles through all characters in turn.
}
\value{
A character vector of the same lenth as \code{x}, but with characters replaced
}
\author{
Andrew Hipp
}
\note{
I have used this function in many projects, and duplicates of it linger in a variety of places.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nuttyName <- "Jimmy-Smith"
notNuttyName <- tidyName("Jimmy-Smith")

## The function is currently defined as
function (x, fixes = c("_", ".", " ", "-", ">"), case = c("lower", 
    "upper", "nochange")) 
{
    x = switch(case[1], lower = tolower(x), upper = toupper(x), 
        nochange = x)
    for (i in fixes) x <- (gsub(i, "", x, fixed = T))
    x
  }
}

