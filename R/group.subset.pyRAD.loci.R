.group.subset.pyRAD.loci <- function(dat, groups, mins = 10, loci = names(dat$DNA), use.tidyName = FALSE, cores = 1, ...) {
# arguments:
#  dat is a subset.pyRAD.loci object
#  groups is a list of individuals (can be named) in the groups
#
# value: a matrix with the number of individuals in each group for each locus

  if(use.tidyName) {
    if(cores > 1) dat$DNA <- mclapply(dat$DNA, function(x) {names(x) <- tidyName(names(x), ...); x}, mc.cores = cores)
	else dat$DNA <- lapply(dat$DNA, function(x) {names(x) <- tidyName(names(x), ...); x})
	groups <- lapply(groups, tidyName, ...)
	}
  if(cores > 1) out <- t(mcmapply(function(x) c(sapply(groups, function(y) sum(names(x) %in% y)), total = length(x)), dat$DNA, mc.cores = cores))
  else out <- t(sapply(dat$DNA, function(x) c(sapply(groups, function(y) sum(names(x) %in% y)), total = length(x))))
  if(!is.na(mins[1])) leave.in <- apply(out, 1, function(x) all(x >= mins))
  else leave.in <- rep(TRUE, dim(out)[1])
  out <- out[leave.in, ]
  if(!is.na(loci[1])) out <- out[intersect(row.names(out), loci), ]
  out
  }
