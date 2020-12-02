**maskPositionsByQuality** - *Mask positions with lo sequencing quality*

Description
--------------------

`maskPositionsByQuality` will replace with an N positions that 
have a sequencing quality score lower that `min_quality`.


Usage
--------------------
```
maskPositionsByQuality(
data,
min_quality = 70,
sequence = "sequence_alignment",
quality_num = "quality_sequence_num"
)
```

Arguments
-------------------

data
:   An AIRR data.frame

min_quality
:   Minimun quality. Positions with sequencing quality 
< `min_qual` will be masked.

sequence
:   Name of the column in `data` with sequence data to be
masked.

quality_num
:   Name of the column in `data` with quality scores (a
string of numeric values, comma separated) that can
be used to mask `sequence`.




Value
-------------------

`data` with one additional field with masked sequences. The 
name of this field is created concatenating `sequence` 
and '_masked'.



Examples
-------------------

```R
### Not run:
maskPositionsByQuality(db, min_quality=90)
```

**Error in checkColumns(data, required_cols)**: object 'db' not found






