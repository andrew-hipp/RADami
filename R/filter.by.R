filter.by <- function(dat, taxa, req,
                      threshInds = 2, threshTaxa = 2,
                      use.tidyName = FALSE, ...) {
  ## returns just loci for which the requested taxa are present at some threshold
  ## default to returning 'all'
  ## updated 2017-01-31 to allow you to subset by a threshold humber of individuals per taxon
  if(!class(dat) %in% c('summary.pyRAD.loci', 'pyRAD.loci')) stop("This function only works with summary.pyRAD.loci datatypes")
  if(class(dat) == 'pyRAD.loci') dat <- dat$radSummary
  if(use.tidyName) {
    row.names(dat$inds.mat) <- tidyName(row.names(dat$inds.mat), ...)
	taxa <- tidyName(taxa, ...)
	}
  dat.mat <- dat$inds.mat[taxa, ] # this is the case if you pass in a summary.pyRAD.loci object
  if(threshold == 'all') threshold <- length(taxa)
  return(names(which(apply(dat.mat, 2, sum) >= threshold)))
  }
