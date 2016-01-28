Alakazam
-------------------------------------------------------------------------------
January 28, 2016  
Version 0.2.8

Lineage, diversity, gene usage and amino acid property analysis R package of 
the Change-O suite.

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