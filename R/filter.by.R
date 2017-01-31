filter.by <- function(dat, taxa, req = NULL,
                      threshInds, threshTaxa = 2,
                      use.tidyName = FALSE, ...) {
  ## returns just loci for which the requested taxa are present at some threshold
  ## default to returning 'all'
  ## updated 2017-01-31 to allow you to subset by a threshold humber of individuals per taxon
  if(!class(dat) %in% c('summary.pyRAD.loci', 'pyRAD.loci')) stop("This function only works with summary.pyRAD.loci datatypes")
  if(class(dat) == 'pyRAD.loci') dat <- dat$radSummary
  if(use.tidyName) {
    row.names(dat$inds.mat) <- tidyName(row.names(dat$inds.mat), ...)
	taxa <- sapply(taxa, tidyName, ...)
	}
  dat.mat <- dat$inds.mat[unlist(unique(taxa)), ]
  tax.thresh.mat <- t(sapply(taxa, function(x) colSums(dat.mat[x, ]) >= threshInds[x]))
  if(length(req) > 0) dat.mat <- dat.mat[, which(colSums(dat.mat[req, ]) == length(req))]
  out <- names(which(colSums(dat.mat) >= threshTaxa))
  return(out)
  }
