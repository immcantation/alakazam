Alakazam
-------------------------------------------------------------------------------
December 18, 2015  
Version 0.2.1

Lineage, diversity, gene usage and amino acid property R package of the 
Change-O suite.

Dependencies
-------------------------------------------------------------------------------
R 3.0  
R packages

  - dplyr
  - ggplot2
  - igraph
  - lazyeval
  - scales
  - seqinr
  - stringi

Mercurial Configuration
-------------------------------------------------------------------------------
Update Mercurial .hgignore file with:  
```
syntax: glob
  .*
  *.Rproj
  man/*.Rd
  inst/doc/*
```

Build Instructions
-------------------------------------------------------------------------------
Install build dependencies:
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

Building with Rstudio:

-  Build -> Configure Build Tools
-  Check use devtools option
-  Check use roxygen option
-  Select configure roxygen options and check everything.
-  Build -> Build and Reload

Building from the R console:

```R
library(roxygen2)
library(devtools)
install_deps()
document()
build(vignettes=FALSE)
install()
```