LG.plot <- function(lociBlast, markerPositions, max.evalue = NA, min.alignL = NA, 
                    lg.name = 'LG', pos.name = 'consensus', lg = NA, 
					tickBounds = c(-0.1, 0.1), label = TRUE, tick.cex = 0.4, text.cex = 0.5, 
					markDupes = c('left', 'right', 'n'), markOverlaps = TRUE, totals = TRUE, 
					return.mat = FALSE,
					...) {
## plots the linkage group locations of loci, given a b6 output of loci blasted to mapped markers, and map positions of the markers
## doesn't really belong in pyRAD, but it's the most convenient spot for it to live right now
## ARGUMENTS:
##  lociBlast - BLASTN of loci against mapped markers, b6 format (http://drive5.com/usearch/manual/blast6out.html)
##  markerPositions - a data.frame with map positions for the targets of the loci blast, with row.names being the marker names (targest of lociBlast)
##  lg.name - column name for the linkage groups indicated in markerPositions
##  pos.name - column name for the map position indicated in markerPostions
##  lg - vector of linkage group names to map

  names(lociBlast) <- c("query", "target", "idPercent", "alignL", "mismatchN", "gapN", 
                        "startPosQuery", "endPosQuery", "startPosTarget", "endPosTarget", 
                        "evalue", "bitscore")
  if(!is.na(max.evalue)) lociBlast <- lociBlast[as.numeric(lociBlast$evalue) <= max.evalue, ]
  if(!is.na(min.alignL)) lociBlast <- lociBlast[as.numeric(lociBlast$alignL) >= min.alignL, ]
  lociBlast <- lociBlast[lociBlast$target %in% row.names(markerPositions), ]
  x <- as.data.frame(cbind(lociBlast, markerPositions[lociBlast$target, ]))
  if(is.na(lg[1])) lg <- sort(unique(x[[lg.name]]))
  lg.ranges <- t(sapply(lg, function(z) range(x[[pos.name]][(x[[lg.name]] == z)], na.rm = TRUE)))
  plot(1, xlim = c(0, length(lg) + 1), ylim = c(-1, max(lg.ranges) + 20), type = 'n', xaxt = 'n', ylab = 'Map distance (cM)', xlab = 'Linkage group', ...)
  axis(1, at = seq(length(lg)), labels = lg, cex.axis = 0.6)
  for(i in 1:length(lg)) {
	segments(i, 0, i, lg.ranges[i, 2], ...)
	x.temp <- x[x[[lg.name]] == lg[i], ]
	segments(i + tickBounds[1], x.temp[[pos.name]], i+tickBounds[2])
	if(label) text(i + tickBounds[2], x.temp[[pos.name]], x.temp$query, pos = 4, cex = tick.cex)
	if(markDupes[1] != 'n') {
	  xpos = switch(markDupes[1], left = i + tickBounds[1] - 0.05, right = i + tickBounds[2] + 0.05)
	  x.dupes <- names(duplicated.mapped.loci(x))
	  message(paste("DUPLICATED LOCI ON", lg[i]))
	  for(j in 1:length(x.dupes)) {
	    ypos <- x.temp[[pos.name]][x.temp$query == x.dupes[j]]
		points(rep(xpos, length(ypos)), ypos, col = j)
        if(sum(x.temp$query == x.dupes[j]) > 0) message(paste("   ", x.dupes[j], "--", sum(x.temp$query == x.dupes[j]), "copies [", paste(x.temp[[pos.name]][x.temp$query == x.dupes[j]], collapse = ', '), "]"))
		} # close for
	  } # close if
	if(markOverlaps) {
	  # browser()
	  x.overlaps <- x.temp[[pos.name]][which(duplicated(x.temp[[pos.name]]))]
	  for(j in unique(x.overlaps)) text(x = i + tickBounds[1], y = j, labels = sum(x.temp[[pos.name]] == j), cex = tick.cex, pos = 2, offset = 0.1)
	  } # close if
	if(totals) {
	  total.queries <- length(unique(x.temp$query))
	  total.targets <- length(unique(x.temp$target))
	  total.positions <- length(unique(x.temp[[pos.name]]))
	  text(i, max(lg.ranges) + c(20,15,10), c(total.queries, total.targets, total.positions), cex = text.cex)
	  } # close if
    } # close i
	if(totals) {
	  text(x = (length(lg) + 0.2), y = (max(lg.ranges) + c(20,15,10)), as.character(c(length(unique(x$query)), length(unique(x$target)), length(unique(x[[pos.name]])))), cex = text.cex, pos = 4)
	  text(0.8, max(lg.ranges) + c(20,15,10), c('RAD loci', 'Contigs', 'Map positions'), cex = text.cex, pos = 2)
	  message(paste("Total queries =", length(unique(x$query))))
	  message(paste("Total targets =", length(unique(x$target))))
	  }
	if(return.mat) return(x)
  } # done
  
duplicated.mapped.loci <- function(x, incomparables = FALSE, ...){
  if(!'mapped.loci' %in% class(x)) warning('We were expecting a mapped.loci object from the matchEm function!')
  dup.q <- unique(x$query[duplicated(x$query)])
  dup.q.list <- lapply(dup.q, function(z) x[x$query == z, ])
  names(dup.q.list) <- dup.q
  dup.q.list
  }