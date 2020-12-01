**getPositionQuality** - *Get a data.frame with sequencing qualities per position*

Description
--------------------

Given a data.frame with sequence quality scores in the form of 
a strings of comma separated numeric values, `getPositionQuality` 
will split the values by ',' and return a data.frame with the values
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
:   An AIRR data.frame

sequence_id
:   Name of the column in `data` with sequence identifiers.

sequence
:   Name of the column in `data` with sequence data.

quality_num
:   Name of the column in `data` with quality scores (as
strings of numeric values, comma separated) for `sequence`.




Value
-------------------

`data` with one additional field with masked sequences. The 
name of this field is created concatenating `sequence` 
and '_masked'.









