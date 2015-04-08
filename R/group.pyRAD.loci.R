group.pyRAD.loci <- function(dat, groups, mins = 10, loci = dimnames(dat$radSummary$inds.mat)[[2]], use.tidyName = TRUE, cores = 8, ...) {
# arguments:
#  dat is a pyRAD.loci object
#  groups is a list of individuals (can be named) in the groups
#  mins is the minimum number of individuals required per group to report on a locus
# value: a matrix with the number of individuals in each group for each locus

  if(use.tidyName) {
    row.names(dat$radSummary$inds.mat) <- tidyName(row.names(dat$radSummary$inds.mat), ...)
	groups <- lapply(groups, tidyName, ...)
	}
  inds.mats <- mclapply(groups, function(x) temp <- dat$radSummary$inds.mat[x[x %in% row.names(dat$radSummary$inds.mat)], ], mc.cores = cores)
  out <- mcmapply(function(x) apply(x, 2, sum), inds.mats, mc.cores = cores)
  if(!is.na(mins[1])) leave.in <- apply(out, 1, function(x) all(x >= mins))
  else leave.in <- rep(TRUE, dim(out)[1])
  out <- cbind(out[leave.in, ], total = colSums(dat$radSummary$inds.mat)[leave.in])
  if(!is.na(loci[1])) out <- out[intersect(row.names(out), loci), ]
  attr(out, 'groups') = groups
  out
  }
 