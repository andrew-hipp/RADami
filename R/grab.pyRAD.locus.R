grab.pyRAD.locus <- function(pyDat, locName, dat.format = c('text', 'fasta', 'matrix')) {
  seqs <- pyDat$seqs[pyDat$locus.index == locName]
  tips <- pyDat$tips[pyDat$locus.index == locName]
  if(dat.format[1] == 'text') names(seqs) <- pyDat$tips[pyDat$locus.index == locName]
  if(dat.format[1] == 'fasta') seqs <- as.character(matrix(c(tips, seqs), 2, byrow = T))
  if(dat.format[1] == 'matrix') seqs <- as.matrix(seqs)
  seqs
}
