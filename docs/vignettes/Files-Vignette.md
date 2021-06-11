# File input and output

As part of the Immcantation suite of tools, the `alakazam` package includes a set of 
built-in functions capable of reading and writing tab-delimited database files created by 
[Change-O](https://changeo.readthedocs.io/en/stable/) into R data.frames. However, due to 
differences in how certain values and sequences are handled, `alakazam::readChangeoDb` and 
`alakazam::writeChangeoDb` will not properly read in [AIRR](https://docs.airr-community.org) 
formatted files. These files should instead be loaded using the functions included 
in the `airr` package (`airr::read_rearrangement` and `airr::write_rearrangement`).

You can read more about how we use both data standards
[here](https://immcantation.readthedocs.io/en/stable/datastandards.html) and 
[here](https://changeo.readthedocs.io/en/stable/standard.html). *Please note that the default 
file format for all functions in Immcantation is the AIRR-C format as of Immcantation 
v4.0.0, which corresponds to alakazam v1.0.0.*

## Reading data

Small example databases for both the Change-O format (`ExampleDbChangeo`) and the AIRR format (`ExampleDb`) 
are included in the `alakazam` package. For specific details about the latter, visit the 
[AIRR Community documentation site](https://docs.airr-community.org/en/stable/datarep/rearrangements.html).


```r
# Set the file paths from inside the package directory
# These files are smaller versions of the example databases previously mentioned
changeo_file <- system.file("extdata", "example_changeo.tab.gz", package="alakazam")
airr_file <- system.file("extdata", "example_airr.tsv.gz", package="alakazam")

# Read in the data
db_changeo <- alakazam::readChangeoDb(changeo_file)
db_airr <- airr::read_rearrangement(airr_file)
```

## Writing data


```r
# Write the data to a tab-delimited file
alakazam::writeChangeoDb(db_changeo, "changeo.tsv")
airr::write_rearrangement(db_airr, "airr.tsv")
```
