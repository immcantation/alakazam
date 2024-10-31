**getPositionQuality** - *Get a data.frame with sequencing qualities per position*

Description
--------------------

`getPositionQuality` takes a data.frame with sequence quality scores 
in the form of a strings of comma separated numeric values, split the quality 
scores values by `","`,  and returns a data.frame with the values
for each position.


Usage
--------------------
```
getPositionQuality(
data,
sequence_id = "sequence_id",
sequence = "sequence_alignment",
quality_num = "quality_alignment_num"
)
```

Arguments
-------------------

data
:   `data.frame` containing sequence data.

sequence_id
:   column in `data` with sequence identifiers.

sequence
:   column in `data` with sequence data.

quality_num
:   column in `data` with quality scores (as
strings of numeric values, comma separated) for `sequence`.




Value
-------------------

`data` with one additional field with masked sequences. The 
name of this field is created concatenating `sequence` 
and '_masked'.



Examples
-------------------

```R
db <- airr::read_rearrangement(system.file("extdata", "test_seq.tsv", package="alakazam"))

```

**Error**: '' does not exist in current working directory ('/home/edelaron/git/Immcantation/alakazam').
```R
fastq_file <- system.file("extdata", "test_seq.fastq", package="alakazam")
db <- readFastqDb(db, fastq_file, quality_offset=-33)

```

**Error in eval(expr, envir, enclos)**: object 'db' not found
```R
head(getPositionQuality(db))
```

**Error in eval(expr, envir, enclos)**: object 'db' not found

See also
-------------------

[readFastqDb](readFastqDb.md) and [maskPositionsByQuality](maskPositionsByQuality.md)






