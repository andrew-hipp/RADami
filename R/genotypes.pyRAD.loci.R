genotypes.pyRAD.loci <- function(dat, groups, loci = 'all', taxa = 'all',
	                             useSnps = c('first', 'all'), concat = c(FALSE, TRUE),
					use.tidyName = FALSE, na.rm = c('none', 'columns', 'rows'), maxAlleles = 2,
					tidyVals = c('-', '.','>', '_', ' '), sortByGroups = TRUE,
					variable.only = FALSE, make.dummy.column = TRUE, alleleDigits = 1, toInteger = TRUE, missingData = "00",
					cores = 1) {
  if(!'subset.pyRAD.loci' %in% class(dat)) stop('Currently, this function is written to require DNAStringSet output from subset.pyRAD.loci,\n
                                                 with only SNPs exported')
  if(taxa[1] != 'all') {
    if(cores > 1) dat$DNA <- mclapply(dat$DNA, function(x) x[names(x) %in% taxa], mc.cores = cores)
    else dat$DNA <- lapply(dat$DNA, function(x) x[names(x) %in% taxa])
    dat$DNA <- dat$DNA[sapply(dat$DNA, length) > 0]
	}
  if(loci[1] != 'all') dat$DNA <- dat$DNA[loci]
  out <- structure(vector('list', length(dat$DNA)), names = names(dat$DNA))
  duplicated.members <- unlist(groups)[duplicated(unlist(groups))]
  if(length(duplicated.members) > 0) warning('Some individuals are duplicated between groups; excluding duplicates from export, including first')
  groups.vector <- structure(integer(length(unlist(groups)[!duplicated(unlist(groups))])), names = unlist(groups)[!duplicated(unlist(groups))])
  for(i in 1:length(groups)) groups.vector[groups[[i]]][!names(groups.vector[groups[[i]]]) %in% duplicated.members] <- i

## 3. Translate SNPs to genotypes
  do.this <- function(y) {
		y.mat <- as.matrix(y)
    if(variable.only) {
	  y.mat.uns <- apply(y.mat, 2, function(x) length(unique(x)) > 1)
      y.mat <- matrix(y.mat[, y.mat.uns], nrow = dim(y.mat)[1], dimnames = list(row.names(y.mat), NULL))
	  }
    if(dim(y.mat)[2] == 0) return('failed')
    if(alleleDigits > 1) pad <- paste(rep("0", alleleDigits-1), collapse = '')
    else pad <- ''
    tip.names <- row.names(y.mat)
    trans.dna <- t(apply(y.mat, 1, function(x) IUPAC_CODE_MAP[x]))
		trans.dna <- t(apply(trans.dna, 1, function(x) gsub('A', paste(pad, '1', sep = ''), x)))
		trans.dna <- t(apply(trans.dna, 1, function(x) gsub('C', paste(pad, '2', sep = ''), x)))
		trans.dna <- t(apply(trans.dna, 1, function(x) gsub('G', paste(pad, '3', sep = ''), x)))
		trans.dna <- t(apply(trans.dna, 1, function(x) gsub('T', paste(pad, '4', sep = ''), x)))
		trans.dna <- t(apply(trans.dna, 1, function(x) {x[nchar(x) == alleleDigits] <- paste(x[nchar(x) == alleleDigits], x[nchar(x) == alleleDigits], sep = ''); return(x)}))
		trans.dna <- as.matrix(trans.dna)[!apply(as.matrix(trans.dna), 1, function(x) any(nchar(x) > maxAlleles * alleleDigits)), ]
		if(na.rm[1] == 'rows') trans.dna <- as.matrix(trans.dna)[!apply(as.matrix(trans.dna), 1, function(x) any(is.na(x))), ]
		if(na.rm[1] == 'columns') trans.dna <- as.matrix(trans.dna)[, !apply(as.matrix(trans.dna), 2, function(x) any(is.na(x)))]
		if(length(trans.dna) == 0) return(0)
		groupMembership <- groups.vector[match(tidyName(row.names(as.matrix(trans.dna)), tidyVals), tidyName(names(groups.vector), tidyVals))]
		if(toInteger) {
	          if(is.null(dim(trans.dna))) trans.dna <- as.integer(trans.dna)
	          else trans.dna <- t(apply(trans.dna, 1, as.integer))
	          }
	        # if(!is.null(dim(trans.dna))) trans.dna <- t(trans.dna)
	        dna.out <- as.data.frame(cbind(groupMembership = groupMembership, trans.dna), stringsAsFactors = FALSE)
		row.names(dna.out) <- tip.names
		if(sortByGroups) dna.out <- dna.out[order(dna.out$groupMembership), ]
		if(make.dummy.column & dim(dna.out)[2] == 2) dna.out$dummy.locus <- rep(11, dim(dna.out)[2])
		return(dna.out)
	}
  out <- mclapply(dat$DNA, function(x) try(do.this(x), silent = TRUE), mc.cores = cores) ## need to figure out what is causing do.this to fail, then get rid of try!
  out <- out[sapply(out, class) %in% c('data.frame', 'matrix')] # gets rid of anything that isn't a matrix or a data.frame
  out <- out[!apply(t(sapply(out, dim)), 1, function(x) sum(x == 0) > 0)] # gets rid of all the matrices in which some dimension == 0
  if(useSnps[1] == 'first') out <- lapply(out, function(x) x[, 1:2])
  groupMembership <- t(sapply(out, function(w) sapply(1:2, function(x) sum(w$groupMembership == x))))
  dimnames(groupMembership)[[2]] <- names(groups)
  if(concat[1]) {
    taxa <- unique(unlist(sapply(out, row.names)))
    out <- do.call(cbind, lapply(out, function(x) x[taxa, 2:(dim(x)[2])]))
    row.names(out) <- taxa
    }
  if(!is.na(missingData)) out[is.na(out)] <- missingData
  attr(out, 'groupMembership') <- groupMembership
  out
}
