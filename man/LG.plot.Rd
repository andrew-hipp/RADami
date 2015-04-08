\name{LG.plot}
\alias{LG.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LG.plot(lociBlast, markerPositions, max.evalue = NA, min.alignL = NA, lg.name = "LG", pos.name = "consensus", lg = NA, tickBounds = c(-0.1, 0.1), label = TRUE, tick.cex = 0.4, text.cex = 0.5, markDupes = c("left", "right", "n"), markOverlaps = TRUE, totals = TRUE, return.mat = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lociBlast}{
%%     ~~Describe \code{lociBlast} here~~
}
  \item{markerPositions}{
%%     ~~Describe \code{markerPositions} here~~
}
  \item{max.evalue}{
%%     ~~Describe \code{max.evalue} here~~
}
  \item{min.alignL}{
%%     ~~Describe \code{min.alignL} here~~
}
  \item{lg.name}{
%%     ~~Describe \code{lg.name} here~~
}
  \item{pos.name}{
%%     ~~Describe \code{pos.name} here~~
}
  \item{lg}{
%%     ~~Describe \code{lg} here~~
}
  \item{tickBounds}{
%%     ~~Describe \code{tickBounds} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{tick.cex}{
%%     ~~Describe \code{tick.cex} here~~
}
  \item{text.cex}{
%%     ~~Describe \code{text.cex} here~~
}
  \item{markDupes}{
%%     ~~Describe \code{markDupes} here~~
}
  \item{markOverlaps}{
%%     ~~Describe \code{markOverlaps} here~~
}
  \item{totals}{
%%     ~~Describe \code{totals} here~~
}
  \item{return.mat}{
%%     ~~Describe \code{return.mat} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (lociBlast, markerPositions, max.evalue = NA, min.alignL = NA, 
    lg.name = "LG", pos.name = "consensus", lg = NA, tickBounds = c(-0.1, 
        0.1), label = TRUE, tick.cex = 0.4, text.cex = 0.5, markDupes = c("left", 
        "right", "n"), markOverlaps = TRUE, totals = TRUE, return.mat = FALSE, 
    ...) 
{
    names(lociBlast) <- c("query", "target", "idPercent", "alignL", 
        "mismatchN", "gapN", "startPosQuery", "endPosQuery", 
        "startPosTarget", "endPosTarget", "evalue", "bitscore")
    if (!is.na(max.evalue)) 
        lociBlast <- lociBlast[as.numeric(lociBlast$evalue) <= 
            max.evalue, ]
    if (!is.na(min.alignL)) 
        lociBlast <- lociBlast[as.numeric(lociBlast$alignL) >= 
            min.alignL, ]
    lociBlast <- lociBlast[lociBlast$target \%in\% row.names(markerPositions), 
        ]
    x <- as.data.frame(cbind(lociBlast, markerPositions[lociBlast$target, 
        ]))
    if (is.na(lg[1])) 
        lg <- sort(unique(x[[lg.name]]))
    lg.ranges <- t(sapply(lg, function(z) range(x[[pos.name]][(x[[lg.name]] == 
        z)], na.rm = TRUE)))
    plot(1, xlim = c(0, length(lg) + 1), ylim = c(-1, max(lg.ranges) + 
        20), type = "n", xaxt = "n", ylab = "Map distance (cM)", 
        xlab = "Linkage group", ...)
    axis(1, at = seq(length(lg)), labels = lg, cex.axis = 0.6)
    for (i in 1:length(lg)) {
        segments(i, 0, i, lg.ranges[i, 2], ...)
        x.temp <- x[x[[lg.name]] == lg[i], ]
        segments(i + tickBounds[1], x.temp[[pos.name]], i + tickBounds[2])
        if (label) 
            text(i + tickBounds[2], x.temp[[pos.name]], x.temp$query, 
                pos = 4, cex = tick.cex)
        if (markDupes[1] != "n") {
            xpos = switch(markDupes[1], left = i + tickBounds[1] - 
                0.05, right = i + tickBounds[2] + 0.05)
            x.dupes <- names(duplicated.mapped.loci(x))
            message(paste("DUPLICATED LOCI ON", lg[i]))
            for (j in 1:length(x.dupes)) {
                ypos <- x.temp[[pos.name]][x.temp$query == x.dupes[j]]
                points(rep(xpos, length(ypos)), ypos, col = j)
                if (sum(x.temp$query == x.dupes[j]) > 0) 
                  message(paste("   ", x.dupes[j], "--", sum(x.temp$query == 
                    x.dupes[j]), "copies [", paste(x.temp[[pos.name]][x.temp$query == 
                    x.dupes[j]], collapse = ", "), "]"))
            }
        }
        if (markOverlaps) {
            x.overlaps <- x.temp[[pos.name]][which(duplicated(x.temp[[pos.name]]))]
            for (j in unique(x.overlaps)) text(x = i + tickBounds[1], 
                y = j, labels = sum(x.temp[[pos.name]] == j), 
                cex = tick.cex, pos = 2, offset = 0.1)
        }
        if (totals) {
            total.queries <- length(unique(x.temp$query))
            total.targets <- length(unique(x.temp$target))
            total.positions <- length(unique(x.temp[[pos.name]]))
            text(i, max(lg.ranges) + c(20, 15, 10), c(total.queries, 
                total.targets, total.positions), cex = text.cex)
        }
    }
    if (totals) {
        text(x = (length(lg) + 0.2), y = (max(lg.ranges) + c(20, 
            15, 10)), as.character(c(length(unique(x$query)), 
            length(unique(x$target)), length(unique(x[[pos.name]])))), 
            cex = text.cex, pos = 4)
        text(0.8, max(lg.ranges) + c(20, 15, 10), c("RAD loci", 
            "Contigs", "Map positions"), cex = text.cex, pos = 2)
        message(paste("Total queries =", length(unique(x$query))))
        message(paste("Total targets =", length(unique(x$target))))
    }
    if (return.mat) 
        return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
