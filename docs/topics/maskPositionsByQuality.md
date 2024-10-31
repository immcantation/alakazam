**maskPositionsByQuality** - *Mask sequence positions with low quality*

Description
--------------------

`maskPositionsByQuality` will replace positions that 
have a sequencing quality score lower that `min_quality` with an
`"N"` character.


Usage
--------------------
```
maskPositionsByQuality(
data,
min_quality = 70,
sequence = "sequence_alignment",
quality_num = "quality_alignment_num"
)
```

Arguments
-------------------

data
:   `data.frame` containing sequence data.

min_quality
:   minimum quality score. Positions with sequencing quality 
less than `min_qual` will be masked.

sequence
:   column in `data` with sequence data to be masked.

quality_num
:   column in `data` with quality scores (a
string of numeric values, comma separated) that can
be used to mask `sequence`.




Value
-------------------

Modified `data` data.frame with an additional field containing 
quality masked sequences. The  name of this field is created 
concatenating the `sequence` name and `"_masked"`.



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
maskPositionsByQuality(db, min_quality=90, quality_num="quality_alignment_num")
```

**Error in eval(expr, envir, enclos)**: object 'db' not found

See also
-------------------

[readFastqDb](readFastqDb.md) and [getPositionQuality](getPositionQuality.md)






