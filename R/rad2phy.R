rad2phy <-
function(pyDat, inds = row.names(pyDat), loci = dimnames(pyDat)[[2]], outfile = 'pyMat.out.phy', padding = 50, verbose = FALSE, logfile = 'rad2phy.log') {
## makes a phylip-style data matrix from rad.mat output, limiting by individuals and loci
  if(class(pyDat) != "rad.mat") warning("I'm expecting output from rad.mat")
  outfile.name = outfile
  outfile = file(outfile, "wt")
  open(outfile)
  cat(paste(length(inds), sum(sapply(pyDat[inds[1], loci], nchar)), "\n"), file = outfile) #header: number of individuals, number of bases
  for(i in inds) {
    if(!verbose) message(paste("Writing DNA line for individual", i))
	cat(i, file = outfile)
	cat(paste(rep(" ", padding - nchar(i)), collapse = ""), file = outfile)
	cat(paste(pyDat[i, loci], collapse = ""), file = outfile)
	cat("\n", file = outfile) # endline
	}
  close(outfile)
  if(!is.na(logfile) & logfile != '') logfile = file(logfile, 'wt')
  open(logfile)
  writeLines(timestamp(), con = logfile)
  writeLines(paste("Filename:", outfile.name), con = logfile)
  writeLines("Loci included in phylip file:", con = logfile)
  writeLines(paste("\t", loci, sep = ''), con = logfile)
  close(logfile)
  }
