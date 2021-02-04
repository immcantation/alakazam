**readFastqDb** - *Load sequencing quality scores from a FASTQ file*

Description
--------------------

`readFastqDb` adds the sequencing quality scores to a data.frame
from a FASTQ file. Matching is done by `sequence_id`.


Usage
--------------------
```
readFastqDb(
data,
fastq_file,
quality_offset = -33,
header = c("presto", "asis"),
sequence_id = "sequence_id",
sequence = "sequence",
sequence_alignment = "sequence_alignment",
v_cigar = "v_cigar",
d_cigar = "d_cigar",
j_cigar = "j_cigar",
np1_length = "np1_length",
np2_length = "np2_length",
v_sequence_end = "v_sequence_end",
d_sequence_end = "d_sequence_end",
style = c("num", "ascii", "both"),
quality_sequence = FALSE
)
```

Arguments
-------------------

data
:   `data.frame` containing sequence data.

fastq_file
:   path to the fastq file

quality_offset
:   offset value to be used by ape::read.fastq. It is 
the value to be added to the quality scores 
(the default -33 applies to the Sanger format and 
should work for most recent FASTQ files).

header
:   FASTQ file header format; one of `"presto"` or 
`"asis"`. Use `"presto"` to specify 
that the fastq file headers are using the pRESTO
format and can be parsed to extract 
the `sequence_id`. Use `"asis"` to skip 
any processing and use the sequence names as they are.

sequence_id
:   column in `data` that contains sequence 
identifiers to be matched to sequence identifiers in 
`fastq_file`.

sequence
:   column in `data` that contains sequence data.

sequence_alignment
:   column in `data` that contains IMGT aligned sequence data.

v_cigar
:   column in `data` that contains CIGAR 
strings for the V gene alignments.

d_cigar
:   column in `data` that contains CIGAR 
strings for the D gene alignments.

j_cigar
:   column in `data` that contains CIGAR 
strings for the J gene alignments.

np1_length
:   column in `data` that contains the number
of nucleotides between the V gene and first D gene 
alignments or between the V gene and J gene alignments.

np2_length
:   column in `data` that contains the number
of nucleotides between either the first D gene and J 
gene alignments or the first D gene and second D gene
alignments.

v_sequence_end
:   column in `data` that contains the 
end position of the V gene in `sequence`.

d_sequence_end
:   column in `data` that contains the 
end position of the D gene in `sequence`.

style
:   how the sequencing quality should be returned;
one of `"num"`, `"phred"`, or `"both"`.
Specify `"num"` to store the quality scores as strings of 
comma separated numeric values. Use `"phred"` to have
the function return the scores as Phred (ASCII) scores. 
Use `"both"` to retrieve both.

quality_sequence
:   specify `TRUE` to keep the quality scores for 
`sequence`. If false, only the quality score
for `sequence_alignment` will be added to `data`.




Value
-------------------

Modified `data` with additional fields:

1.  `quality_alignment`:     A character vector with ASCII Phred 
scores for `sequence_alignment`.
1.  `quality_alignment_num`: A character vector, with comma separated 
numerical quality values for each 
position in `sequence_alignment`.
1.  `quality`:      A character vector with ASCII Phred 
scores for `sequence`.
1.  `quality_num`:  A character vector, with comma separated 
numerical quality values for each 
position in `sequence`.




Examples
-------------------

```R
db <- airr::read_rearrangement(system.file("extdata", "test_seq.tsv", package="alakazam"))
fastq_file <- system.file("extdata", "test_seq.fastq", package="alakazam")
db <- readFastqDb(db, fastq_file, quality_offset=-33)
```



See also
-------------------

[maskPositionsByQuality](maskPositionsByQuality.md) and [getPositionQuality](getPositionQuality.md)






