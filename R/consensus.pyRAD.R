consensus.pyRAD <-
function(pyIn, from = NA, to = NA, fastaNames = T, writeFile = 'rads.con.txt', cores = 1, ...) {
## originally used seqinr to generate a consensus sequence for each pyRAD locus
## 2013-01-04: updated to use Biostrings, which works better 
##             - deleted arguments: method = 'majority', threshold = 0.001
  if(class(pyIn) != "pyRAD.loci") stop("pyRAD input required, from read.pyRAD")
  allLoci <- unique(as.character(pyIn$locus.index))
  allLoci <- allLoci[!allLoci == ""] ## this should probably be moved to read.pyRAD
  if(!is.na(from)) allLoci <- allLoci[from:to]
  seqs <- as.character(pyIn$seqs)
  loc.index <- as.character(pyIn$locus.index)
  if(cores == 1) {
    out <- character(0)
    for (i in allLoci) out <- c(out, consensusString(DNAStringSet(gsub("-", "N", seqs[loc.index == i])), ...))
    }
  else out <- mcmapply(function(i) consensusString(DNAStringSet(gsub("-", "N", seqs[loc.index == i])), ...), allLoci, mc.cores = cores)
  if(fastaNames) allLoci <- paste(">", allLoci, sep = "")
  names(out) <- allLoci
  if(!is.na(writeFile)) write.table(out, writeFile, sep = "\n", quote = F, col.names = F)
  return(out)
  }
