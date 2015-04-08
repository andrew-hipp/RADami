tidyName <- function(x, fixes = c('_', '.', ' ', '-', '>'), case = c('lower', 'upper', 'nochange')) {
  x = switch(case[1], lower = tolower(x), upper = toupper(x), nochange = x)
  for(i in fixes) x <- (gsub(i, "", x, fixed = T))
  x
  }
