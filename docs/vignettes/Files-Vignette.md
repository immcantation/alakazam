# File input and output
As part of the Immcantation suite of tools, the `alakazam` package includes a set of 
built-in functions capable of reading and writing tab-delimited database files created by 
[Change-O](https://changeo.readthedocs.io/en/stable/) into R data.frames. However, due to 
differences in how certain values and sequences are handled, `alakazam::readChangeoDb` and 
`alakazam::writeChangeoDb` will not properly read in [AIRR](https://docs.airr-community.org/en/stable) 
formatted files. These files should instead be loaded using the functions included in the `airr` package 
(`airr::read_rearrangement` and `airr::write_rearrangement`).

You can read more about how we use both data standards
[here](https://immcantation.readthedocs.io/en/stable/datastandards.html) and 
[here](https://changeo.readthedocs.io/en/stable/standard.html). *Please note that the default 
file format for all functions in Immcantation is the AIRR-C format as of release 4.0.0.*

## Example data

Small example databases for both the Change-O format (`ExampleDbChangeo`) and the AIRR format (`ExampleDb`) 
are included in the `alakazam` package. For specific details about the latter, visit the 
[AIRR Community documentation site ](https://docs.airr-community.org/en/latest/datarep/rearrangements.html#fields).


```r
# Load required packages
library(airr)
library(alakazam)

# Read in the data
db_changeo <- readChangeoDb("../data-raw/ExampleDbChangeo.gz")
db_airr <- read_rearrangement("../data-raw/ExampleDb.gz")

# Write the data to a tab-delimited file
writeChangeoDb(db_changeo, "changeo.tsv")
write_rearrangement(db_airr, "airr.tsv")
```