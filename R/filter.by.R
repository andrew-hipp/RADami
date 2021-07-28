filter.by <- function(dat, taxa, req = NULL,
                      threshInds = 3, threshTaxa = NULL,
                      use.tidyName = FALSE, ...) {
  ## returns just loci for which the requested taxa are present at some threshold
  ## default to returning 'all'
  ## updated 2017-01-31 to allow you to subset by a threshold humber of individuals per taxon
  if(!is.null(req)) warning('req deprecated; ignored in this version')
  if(!is.null(threshTaxa)) warning('threshTaxa deprecated; ignored in this version')
  if(length(threshInds) > 1) {
    warning('threshInds limited to one value; using threshInds[1]')
    threshInds <- threshInds[1]
  }
  if(!class(dat) %in% c('summary.pyRAD.loci', 'pyRAD.loci')) stop("This function only works with summary.pyRAD.loci datatypes")
  if(class(dat) == 'pyRAD.loci') dat <- dat$radSummary
  if(use.tidyName) {
    row.names(dat$inds.mat) <- tidyName(row.names(dat$inds.mat), ...)
	  taxa <- sapply(taxa, tidyName, ...)
	}
  # if(length(threshInds) == 1) threshInds <- structure(rep(threshInds, length(taxa)), names = taxa)
  dat.mat <- dat$inds.mat[unlist(unique(taxa)), ]
  # tax.thresh.mat <- sapply(taxa, function(x) colSums(dat.mat[x, ]) >= threshInds[x])
  # if(!is.null(req)) {
  #   if(length(req) > 1) loc.temp <- names(which(colSums(tax.thresh.mat[req, ]) == length(req)))
  #   if(length(req) == 1) loc.temp <- names(which(tax.thresh.mat[req, ]))
  # } else loc.temp <- colnames(tax.thresh.mat)
  out <- names(which(colSums(dat.mat[taxa, ]) >= threshInds))
  return(out)
  }

# filter.by(rads, threshInds = 3, taxa = row.names(rads$radSummary$inds.mat))
