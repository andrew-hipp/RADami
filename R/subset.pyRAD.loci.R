subset.pyRAD.loci <-
function(x, loci = colnames(x$radSummary$inds.mat), taxa = row.names(x$radSummary$inds.mat),
         format = 'DNAStringSet', reportInterval = 500, mins = 1,
         nucVarType = c("verystrict", "strict", "relaxed"),
		 use.tidyName = FALSE,
         snpsOnly = FALSE, cores = 1, ...) {
## Arguments:
##  x = pyRAD.loci object
##  loci = loci to use (by name)
##  taxa = taxa to use (by name)
##  format = only DNAStringSet export supported now
##  reportInterval = interval to use for reporting subsetting progress
##  nucVarType = either strict for requiring variability to be due to only unambiguous nucleotides (strict, verystrict) or allowing ambiguities to encode variability (relaxed), and insisting on parsimony informativeness (verystrict) or just variability (strict, relaxed) at at least one site
##  2014-03-18: added snpLocs to out to return where the snps are
  excludedNucs <- switch(nucVarType[1], verystrict = 5:17, strict = 5:17, relaxed = 15:17)
  nucThresh <- switch(nucVarType[1], verystrict = 2, strict = 1, relaxed = 1)
  if(use.tidyName) inds.vector <- tidyName(x$tips, ...) %in% tidyName(taxa, ...)
  else inds.vector <- x$tips %in% taxa

  do.it.dna <- function(i) {
	  seq.index <- x$locus.index == i & inds.vector
	  dna.temp <- DNAStringSet(x$seqs[seq.index])
	  names(dna.temp) <- x$tips[seq.index]
	  return(dna.temp)
	  } # close do.it.dna
  if(cores != 1) out <- list(DNA = mclapply(loci, do.it.dna, mc.cores = cores))
  else out <- list(DNA = lapply(loci, do.it.dna))
  names(out$DNA) <- loci
  if(cores != 1) out$snpLocs <- mclapply(out$DNA, function(i) try(which(apply(consensusMatrix(i)[-c(excludedNucs), ], 2, function(x) sum(x >= nucThresh) > 1))), mc.cores = cores)
  else out$snpLocs <- lapply(out$DNA, function(i) try(which(apply(consensusMatrix(i)[-c(excludedNucs), ], 2, function(x) sum(x >= nucThresh) > 1))))
  if(snpsOnly) {
	do.it.snps <- function(i) {
	    tempDNA <- as.matrix(as.matrix(out$DNA[[i]])[, out$snpLocs[[i]]])
		snps.out <- try(DNAStringSet(apply(tempDNA, 1, paste, collapse = '')))
		if(class(snps.out) == 'try-error') message(paste('failed on locus', i))
		return(snps.out)
	    } # close do.it.snps
    if(cores != 1) out$DNA <- mclapply(loci, do.it.snps, mc.cores = cores)
    else out$DNA <- lapply(loci, do.it.snps)
    names(out$DNA) <- loci
    out$DNA <- out$DNA[sapply(out$DNA, class) != 'try-error']
    } # close if snpsOnly
  the.haves <-  which(sapply(out$DNA, function(x) max(width(x)) > 0))
  # print(the.haves)
  out$DNA <- out$DNA[the.haves]
  out$ntaxa <- sapply(out$DNA, length)
  out$variable <- sapply(out$snpLocs, function(i) length(i) > 0)
  class(out) <- 'subset.pyRAD.loci'
  return(out)
  }
