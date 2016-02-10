Alakazam
-------------------------------------------------------------------------------
February 10, 2016  
Version 0.2.3

Lineage, diversity, gene usage and amino acid property analysis R package of 
the Change-O suite.

Dependencies
-------------------------------------------------------------------------------
R 3.1.2  
R packages

  - dplyr
  - ggplot2
  - igraph
  - lazyeval
  - scales
  - seqinr
  - stringi
  - knitr
  - rmarkdown

Build Instructions
-------------------------------------------------------------------------------
Install build dependencies:
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

Building with Rstudio:

- _Build_ -> _Configure Build Tools_
- Check the _Use devtools package functions_ option
- Check the _Generate documentation with Roxygen_ option
- Select _Configure..._ Roxygen options and check everything.
- _Build_ -> _Build and Reload_

Building from the R console:

```R
devtools::install_deps()
devtools::document()
devtools::build()
devtools::install()
```

Optionally, you can skip the vignettes:
```R
devtools::build(vignettes=FALSE)
```
