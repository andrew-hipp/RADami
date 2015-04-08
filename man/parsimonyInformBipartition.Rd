\name{parsimonyInformBipartition}
\alias{parsimonyInformBipartition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parsimonyInformBipartition(dat, bipartition, minPerGroup = 5, return.option = c("max", "mean", "first", "all", "mean.all"), use.tidyNames = TRUE, remove.errors = FALSE, cores = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{bipartition}{
%%     ~~Describe \code{bipartition} here~~
}
  \item{minPerGroup}{
%%     ~~Describe \code{minPerGroup} here~~
}
  \item{return.option}{
%%     ~~Describe \code{return.option} here~~
}
  \item{use.tidyNames}{
%%     ~~Describe \code{use.tidyNames} here~~
}
  \item{remove.errors}{
%%     ~~Describe \code{remove.errors} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dat, bipartition, minPerGroup = 5, return.option = c("max", 
    "mean", "first", "all", "mean.all"), use.tidyNames = TRUE, 
    remove.errors = FALSE, cores = 2) 
{
    require(plyr)
    if (!"subset.pyRAD.loci" \%in\% class(dat)) 
        stop("this function requires an object of class subset.pyRAD.loci")
    if (use.tidyNames) 
        bipartition = lapply(bipartition, tidyName)
    nucs = c("a", "g", "c", "t", "A", "G", "C", "T")
    mat.stats <- function(datmat, bip) {
        mat <- as.matrix(datmat[row.names(datmat) \%in\% bip, ])
        mat.sums <- lapply(apply(mat, 2, count), function(x) x[x$x \%in\% 
            nucs, ])
        dom.mat <- cbind(do.call(rbind, lapply(mat.sums, function(y) y[which(y$freq == 
            max(y$freq))[1], ])), total = sapply(mat.sums, function(x) sum(x$freq)))
        return(dom.mat)
    }
    do.it <- function(workingMat, option = c("mean", "first", 
        "all", "mean.all")) {
        workingMat <- as.matrix(workingMat)
        N = dim(workingMat)[2]
        variable <- apply(as.character(phyDat(workingMat)), 2, 
            function(x) length(unique(x[x \%in\% nucs]))) > 1
        if (sum(variable) == 0) 
            return(0)
        workingMat <- as.matrix(workingMat[, variable])
        if (use.tidyNames) 
            row.names(workingMat) <- tidyName(row.names(workingMat))
        dom.mat1 <- mat.stats(workingMat, bipartition[[1]])
        dom.mat2 <- mat.stats(workingMat, bipartition[[2]])
        if (any(max(dom.mat1$total) < minPerGroup, max(dom.mat2$total) < 
            minPerGroup)) 
            return(NA)
        statNum <- dom.mat1$freq + dom.mat2$freq
        statNum <- ifelse(as.character(dom.mat1$x) == as.character(dom.mat2$x), 
            dom.mat1$total + dom.mat2$total - statNum, statNum)
        stat <- statNum/(dom.mat1$total + dom.mat2$total)
        if (option[1] == "max") 
            out <- max(stat, na.rm = T)
        if (option[1] == "first") 
            out <- stat[1]
        if (option[1] == "all") 
            out <- stat
        if (option[1] == "mean") 
            out <- mean(stat, na.rm = T)
        if (option[1] == "mean.all") 
            out <- sum(stat, na.rm = T)/N
        return(out)
    }
    out <- mclapply(dat$DNA, do.it, option = return.option[1], 
        mc.cores = cores)
    if (return.option[1] \%in\% c("max", "first", "mean", "mean.all")) 
        out <- unlist(out)
    if (remove.errors) 
        out <- out[-grep("Error", out)]
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
