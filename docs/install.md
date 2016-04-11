Download
-------------------------------------------------------------------------------

The latest stable release of Alakazam can be downloaded from 
[CRAN](http://cran.rstudio.com/web/packages/alakazam) or 
[Bitbucket](https://bitbucket.org/kleinstein/alakazam/downloads).

Installing Released Versions
-------------------------------------------------------------------------------

The simplest way to install Alakazam is via CRAN:

```R
install.packages("alakazam")
```

Downloaded source builds from Bitbucket may be installed in the usual way:

```R
install.packages("alakazam_x.y.z.tar.gz", repos=NULL, type="source")
```

Building Development Versions
-------------------------------------------------------------------------------

To build from the [source code](http://bitbucket.org/kleinstein/alakazam),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

To install the latest development code via devtools:

```R
library(devtools)
install_bitbucket("kleinstein/alakazam@default")
```

Note, using `install_bitbucket` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using devtools, 
roxygen and knitr:

```R
library(devtools)
install_deps()
document()
build()
install()
```
