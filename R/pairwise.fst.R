pairwise.fst.loci <- function(dat, group.list, to.do,
                              minInds = 3, do.hist = FALSE, cores = 1,
                              ...) {
  fst <- gt <- gr <- list()
  sp.pairs <- vector('list', 0)
  for(i in 2:length(to.do)) {
    for(j in seq(i-1)) {
      message("DOING:")
	  sp.pairs[[length(sp.pairs) + 1]] <- to.do[c(i, j)]
	  message(paste(c(to.do[i], to.do[j]), collapse = ' - '))
	  gr[[to.do[i]]][[to.do[j]]] <- .group.subset.pyRAD.loci(dat, group.list[c(to.do[i], to.do[j])],
                                                           mins = minInds, cores = cores)
	  gt[[to.do[i]]][[to.do[j]]] <- genotypes.pyRAD.loci(dat, group.list[c(to.do[i], to.do[j])],
                                                       taxa = unlist(group.list[c(to.do[i], to.do[j])]),
                                                       loci = row.names(gr[[to.do[i]]][[to.do[j]]]),
                                                       cores = cores, na.rm = 'columns', ...)
      fst[[to.do[i]]][[to.do[j]]] <- mclapply(gt[[to.do[i]]][[to.do[j]]],
                                              function(x) try(wc(x), silent = TRUE),
                                              mc.cores=cores)
      fst[[to.do[i]]][[to.do[j]]] <- fst[[to.do[i]]][[to.do[j]]][-which(sapply(fst[[to.do[i]]][[to.do[j]]], function(x) class(x) == 'try-error'))]
	  if(do.hist) {
	    pdf(paste('fst', to.do[i], to.do[j], 'pdf', format(Sys.time(), "%Y-%m-%d-%H.%M"), sep = '.'))
	    hist(sapply(fst[[to.do[i]]][[to.do[j]]], function(x) x$FST), 20,
           xlab = 'FST', main = paste('Pairwise FST,', to.do[i], '-', to.do[j]))
	    dev.off()
		}
	  } #close j
	} #close i
	out <- list(fst = fst, groups = gr, genotypes = gt, sp.pairs = sp.pairs,
              snpLocs = dat$snpLocs, timestamp = Sys.time())
	return(out)
	}
