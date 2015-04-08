compare.loci <-
function(taxa, radMat) {
  ## taxa is a pair of taxon vectors
  ## radMat is an inds.mat
  if(length(taxa[[1]]) > 1) v1 <- apply(radMat[taxa[[1]], ], 2, any)
  else v1 <- radMat[taxa[[1]], ]
  if(length(taxa[[2]]) > 1) v2 <- apply(radMat[taxa[[2]], ], 2, any)
  else v2 <- radMat[taxa[[2]], ]
  if(length(taxa) == 3) {
    if(length(taxa[[3]]) > 1) v3 <- apply(radMat[taxa[[3]], ], 2, any)
    else v3 <- radMat[taxa[[2]], ]
	v.mat <- do.call(rbind, list(v1, v2, v3))
    }
  else v.mat <- rbind(v1, v2)
  v.out <- sum(apply(v.mat, 2, all))
  v.denom <- sum(apply(v.mat, 2, any))
  out <- c(abs = v.out, prop = v.out / v.denom)
  return(out)
}

do.all.nodes <-
function(tr, rads, analysis) {
  to.do <- unique(tr$edge[, 1])
  out <- t(sapply(to.do, function(x) compare.loci(node.tips(tr, x, analysis), rads)))
  row.names(out) <- to.do
  out 
}

node.tips <-
function(tr, node, getOut = c('desc.og', 'desc'), add.og = '>AC_h') {
  to.do <- tr$edge[which(tr$edge[,1] == node), 2]
  out <- lapply(Descendants(tr, to.do, 'tips'), function(x) tr$tip.label[x])
  if(getOut[1]== 'desc.og') out[[3]] <- c(tr$tip.label[!tr$tip.label %in% unique(unlist(out))], add.og)
  return(out)
}

do.a.rad.comparison.tree <- function(tr, rads, cex.scalar = 2, analysis = 'desc.og', ...) {
  nodes.comparison <- do.all.nodes(tr, rads, analysis)
  plot(tr)
  nodelabels(pie = nodes.comparison[,'prop'], node = as.numeric(row.names(nodes.comparison)), cex = cex.scalar * nodes.comparison[, 'abs']/max(nodes.comparison[, 'abs']), piecol = c('black', 'white'), ...)
  return(nodes.comparison)
  }