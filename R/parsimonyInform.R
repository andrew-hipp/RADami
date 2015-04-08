parsimonyInformBipartition <- function(dat, bipartition, minPerGroup = 5, 
     return.option = c('max', 'mean', 'first', 'all', 'mean.all'), use.tidyNames = TRUE, remove.errors = FALSE, cores = 2) {
## calculates the parsimony informativeness of a locus / dataset for one bipartition
## Arguments:
##   dat = an object of class subset.pyRAD.loci
##   bipartition = list of two sets of tips representing the bipartition
## currently conservative -- only uses unambiguous sites
## 

## calculate per site as ((number tips with dominant nucleotide for set 1) 
##                       + (number of tips with dominant nucleotide for set 2)
##                       - (total tips if dominant in set 1 is the same as dominant in set 2))
##                       / total number of tips
require(plyr) #CHANGE THIS TO importfrom(plyr, count)

## per locus, calculate so that it ranges from 0 to 1, where a zero is only when there are no variable sites, 
## and a 1 is when one or all variable nucleotides are perfect. Then, do stats either over all SNPs or just for first SNP in each locus
  if(!"subset.pyRAD.loci" %in% class(dat)) stop("this function requires an object of class subset.pyRAD.loci") 
  if(use.tidyNames) bipartition = lapply(bipartition, tidyName)
  nucs = c('a', 'g', 'c', 't', 'A', 'G', 'C', 'T')
  ## nasty embedded function
  mat.stats <- function(datmat, bip) { 
	mat <- as.matrix(datmat[row.names(datmat) %in% bip, ]) # use as.matrix here to ensure that we don't get a vector, when only one nucleotide is variable
	mat.sums <- lapply(apply(mat, 2, count), function(x) x[x$x %in% nucs, ])
	# in dom.mat below, if there is a tie for most common nucleotide, the first is taken; 
	# it should not matter on average whether the most common nucleotide matches the other group most common in this case, 
	# as 0.5 == 1-0.5
	dom.mat <- cbind(do.call(rbind, lapply(mat.sums, function(y) y[which(y$freq == max(y$freq))[1], ])), total = sapply(mat.sums, function(x) sum(x$freq)))
	return(dom.mat)
	}
  do.it <- function(workingMat, option = c('mean', 'first', 'all', 'mean.all')) {
	workingMat <- as.matrix(workingMat) ## need to pass along a matrix
    N = dim(workingMat)[2]
	variable <- apply(as.character(phyDat(workingMat)),2, function(x) length(unique(x[x %in% nucs]))) > 1
    if(sum(variable) == 0) return(0) # even if we get past this without returning 0, there may be columns that have ambiguities
	workingMat <- as.matrix(workingMat[, variable]) # use as.matrix here to ensure that we don't get a vector, when only one nucleotide is variable
	if(use.tidyNames) row.names(workingMat) <- tidyName(row.names(workingMat))
    dom.mat1 <- mat.stats(workingMat, bipartition[[1]])
	dom.mat2 <- mat.stats(workingMat, bipartition[[2]])
	if(any(max(dom.mat1$total) < minPerGroup, max(dom.mat2$total) < minPerGroup)) return(NA) # return NA if we don't have a minimum number of individuals for at least some columns
	statNum <- dom.mat1$freq + dom.mat2$freq
	statNum <- ifelse(as.character(dom.mat1$x) == as.character(dom.mat2$x), dom.mat1$total + dom.mat2$total - statNum, statNum)
	stat <- statNum / (dom.mat1$total + dom.mat2$total)
	if(option[1] == 'max') out <- max(stat, na.rm = T)
	if(option[1] == 'first') out <- stat[1]
	if(option[1] == 'all') out <- stat
	if(option[1] == 'mean') out <- mean(stat, na.rm = T)
	if(option[1] == 'mean.all') out <- sum(stat, na.rm = T) / N
	return(out)
	}
  out <- mclapply(dat$DNA, do.it, option = return.option[1], mc.cores = cores)
  if(return.option[1] %in% c('max', 'first', 'mean', 'mean.all')) out <- unlist(out)
  if(remove.errors) out <- out[-grep('Error', out)] ## I don't think this is actually needed when the data are formatted correctly
  out
  }