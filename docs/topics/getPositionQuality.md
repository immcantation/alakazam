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
db <- airr::read_rearrangement(system.file("extdata", "example_quality.tsv", package="alakazam"))
fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")
db <- readFastqDb(db, fastq_file, quality_offset=-33)
head(getPositionQuality(db))

```

*Warning*:NAs introduced by coercion
```
  position quality_alignment_num       sequence_id nt
1        1                    90 CGCTTTTCGGATTGGAA  C
2        2                    90 CGCTTTTCGGATTGGAA  A
3        3                    90 CGCTTTTCGGATTGGAA  G
4        4                    90 CGCTTTTCGGATTGGAA  C
5        5                    90 CGCTTTTCGGATTGGAA  T
6        6                    90 CGCTTTTCGGATTGGAA  G

```



See also
-------------------

[readFastqDb](readFastqDb.md) and [maskPositionsByQuality](maskPositionsByQuality.md)






