filter.by <- function(dat, taxa, req = NULL,
                      threshInds = 3, threshTaxa = 2,
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
  if(length(threshInds) == 1) threshInds <- structure(rep(threshInds, length(taxa)), names = names(taxa))
  dat.mat <- dat$inds.mat[unlist(unique(taxa)), ]
  tax.thresh.mat <- t(sapply(names(taxa), function(x) colSums(dat.mat[taxa[[x]], ]) >= threshInds[x]))
  if(length(req) > 1) loc.temp <- names(which(colSums(tax.thresh.mat[req, ]) == length(req)))
  else if (length(req) == 1) loc.temp <- names(which(tax.thresh.mat[req, ]))
  else loc.temp <- colnames(tax.thresh.mat)
  out <- list(loci = intersect(loc.temp, names(which(colSums(tax.thresh.mat) >= threshTaxa))),
              loc.mat = tax.thresh.mat)
  return(out)
  }
