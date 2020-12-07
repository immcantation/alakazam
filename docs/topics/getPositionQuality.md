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









