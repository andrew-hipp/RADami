write.DNAStringSet <-
function(x, format= c('phylip', 'fasta'), padding = 30, filename = 'DNAStringSetOut.phy', fastaPrefix = '>') {
  # writes a sequence matrix to phylip or fasta format
  x.width <- width(x)[1]
  x <- as.character(x)
  if(format[1] == 'phylip') {
    for(i in 1:length(x)) x[i] <- paste(names(x)[i], paste(rep(" ", (padding - nchar(names(x)[i]))), collapse = ''), x[i], sep = '')
    writeLines(c(paste(length(x), x.width), x), filename)
	}
  if(format[1] == 'fasta') {
    out <- as.character(matrix(c(paste(fastaPrefix, names(x), sep = ''), x), nrow = 2, byrow = T))
	writeLines(out, filename)
	}
  return(0)
  }
