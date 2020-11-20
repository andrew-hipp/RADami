# RADami
RADseq helper R package

To install, first get the IRanges and BioStrings packages:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IRanges")
BiocManager::install("Biostrings")
```

Then, install using devtools:

```
require('devtools')
install_github('andrew-hipp/RADami')
```

Installation has been tested on Windows 10 in R 4.0.3, "Bunny-Wunnies Freak Out" (fails in R 4.0.2)
