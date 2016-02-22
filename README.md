Alakazam - Immunoglobulin clonal lineage and diversity analysis
-------------------------------------------------------------------------------

The alakazam R package provides a set of tool to investigate lymphocyte receptor 
clonal lineages, diversity and other repertoire level properties, with a focus on 
high-throughput immunoglobulin (Ig) sequencing.

Alakazam serves five main purposes:

1. Providing core functionality for other R packages in the Change-O suite. This
   includes common tasks such as file I/O, basic DNA sequence manipulation, and
   interacting with V(D)J segment and gene annotations.
2. Providing an R interface for interacting with the output of the pRESTO 
   tool suite.
3. Performing lineage reconstruction on clonal populations of immunoglobulin 
   (Ig) sequences. 
4. Performing clonal abundance and diversity analysis on lymphocyte repertoires.
5. Performing physicochemical property analyses of lymphocyte receptor sequences.

Build Instructions
-------------------------------------------------------------------------------
Install build dependencies:
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

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
