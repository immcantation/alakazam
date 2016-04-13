





**readChangeoDb** - *Read a Change-O tab-delimited database file*

Description
--------------------

`readChangeoDb` reads a tab-delimited database file created by a Change-O tool 
into a data.frame.


Usage
--------------------
```
readChangeoDb(file, select = NULL, drop = NULL, seq_upper = TRUE)
```

Arguments
-------------------

file
:   tab-delimited database file output by a Change-O tool.

select
:   columns to select from database file.

drop
:   columns to drop from database file.

seq_upper
:   if `TRUE` convert sequence columns to upper case;
if `FALSE` do not alter sequence columns. See Value 
for a list of which columns are effected.



Value
-------------------

A data.frame of the database file. Columns will be imported as is, except for 
the following columns which will be explicitly converted into character 
values:

+ SEQUENCE_ID
+ CLONE
+ SAMPLE

And the following sequence columns which will converted to upper case if
`seq_upper=TRUE` (default).

+ SEQUENCE_INPUT
+ SEQUENCE_VDJ
+ SEQUENCE_IMGT
+ JUNCTION
+ GERMLINE_IMGT
+ GERMLINE_IMGT_D_MASK




Examples
-------------------

```R
# Load example data
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")

# Read all columns in and convert sequence fields to upper case
df <- readChangeoDb(file)

# Subset columns and convert sequence fields to upper case
df <- readChangeoDb(file, select=c("SEQUENCE_ID", "SEQUENCE_IMGT"))

# Drop columns and do not alter sequence field case
df <- readChangeoDb(file, drop=c("D_CALL", "DUPCOUNT"), seq_upper=FALSE)
```



See also
-------------------

Wraps [read.table](http://www.inside-r.org/r-doc/utils/read.table) and [fread](http://www.inside-r.org/packages/cran/data.table/docs/fread). 
See [writeChangeoDb](writeChangeoDb.md) for writing to Change-O files.



