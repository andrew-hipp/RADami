rad2nex <-
  function(pyDat, inds = row.names(pyDat), indNames = NA,
          fillBlanks = NA, fillChar = '-',
         loci = dimnames(pyDat)[[2]], outfile = 'pyMat.out.nex',
         verbose = FALSE, logfile = 'rad2nex.log', ...) {
## makes a nexus-style data matrix from rad.mat output,
##   limiting by individuals and loci
  if(class(pyDat) != "rad.mat") warning("I'm expecting output from rad.mat")
  temp <- apply(pyDat[inds, loci], 1, paste, collapse = '')
  if(!is.na(indNames[1])) names(temp) <- indNames
  if(!is.na(fillBlanks[1]))
    temp <- c(temp,
      structure(rep(paste(rep(fillChar, nchar(temp[1])),
                              collapse = ''),
                    length(fillBlanks)),
                names = fillBlanks)
              ) # close rbind
  if(verbose) {
    message("Writing nexus file")
    message(paste("Substitute this number for nchar in your file:", nchar(temp[1])))
  }
  write.nexus.data(temp, outfile, ...)
  if(!is.na(logfile) & logfile != '') logfile = file(logfile, 'wt')
  open(logfile)
  writeLines(timestamp(), con = logfile)
  writeLines(paste("Filename:", outfile), con = logfile)
  writeLines(paste("Number of characters:", nchar(temp[1])), con = logfile)
  writeLines("Loci included in nexus file:", con = logfile)
  writeLines(paste("\t", loci, sep = ''), con = logfile)
  close(logfile)
  }
