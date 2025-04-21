# Download and installation

Download
-------------------------------------------------------------------------------

The latest stable release of `alakazam` can be downloaded from 
<a href="http://cran.rstudio.com/web/packages/alakazam" target="_blank">CRAN</a>
or <a href="https://github.com/immcantation/alakazam/tags" target="_blank">GitHub</a>.

Installing Released Versions
-------------------------------------------------------------------------------

The simplest way to install `alakazam` is via CRAN:

```R
install.packages("alakazam")
```

Downloaded source builds from GitHub may be installed in the usual way:

```R
install.packages("alakazam_x.y.z.tar.gz", repos = NULL, type = "source")
```

If you have any trouble installing the package, it may be due to the Bioconductor 
dependencies. You can run the following command to see what other packages may be needed:

```R
available.packages()["alakazam", "Imports"]
```

Alternatively, you can use Bioconductor's `install` function:

```R
install.packages("BiocManager")
BiocManager::install("alakazam")
```

Building Development Versions
-------------------------------------------------------------------------------

To build from the [source code](http://github.com/immcantation/alakazam),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools:

```R
library(devtools)
install_github("immcantation/alakazam@master")
```

Note, using `install_github` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using `devtools`, 
`roxygen` and `knitr`:

```R
library(devtools)
install_deps()
document()
build()
install()
```

Some users might experience issues with alakazam if certain dependencies are not installed correctly. One straightforward way to ensure all required packages are available is to install alakazam via Bioconductor:

```{r}
install.packages("BiocManager")
BiocManager::install("alakazam")
```

If you still encounter issues with missing dependencies, you can use the [Posit Public Package Manager](https://packagemanager.posit.co/client/#/) to identify and install them manually. 


If problems persist while installing alakazam, feel free to contact us [here](https://immcantation.readthedocs.io/en/stable/about.html) for support.
