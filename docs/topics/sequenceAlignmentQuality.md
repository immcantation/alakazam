**sequenceAlignmentQuality** - *Retrieve sequencing quality scores from tabular data*

Description
--------------------

`sequenceAlignmentQuality` is used internally by `readFastqDb` to 
process the sequencing quality scores loaded from a `fastq` file.


Usage
--------------------
```
sequenceAlignmentQuality(
data,
sequence_id = "sequence_id",
sequence = "sequence",
sequence_alignment = "sequence_alignment",
quality_sequence = "quality_sequence",
quality_sequence_num = "quality_sequence_num",
v_cigar = "v_cigar",
d_cigar = "d_cigar",
j_cigar = "j_cigar",
np1_length = "np1_length",
np2_length = "np2_length",
v_sequence_end = "v_sequence_end",
d_sequence_end = "d_sequence_end",
raw = FALSE
)
```

Arguments
-------------------

data
:   `data.frame` containing sequence data.

sequence_id
:   column in `data` that contains sequence 
identifiers to be matched to sequence identifiers in 
`fastq_file`.

sequence
:   column in `data` that contains sequence data.

sequence_alignment
:   column in `data` that contains 
IMGT aligned sequence data.

quality_sequence
:   column in `data` that contains 
sequencing quality as Phred scores.

quality_sequence_num
:   column in `data` that contains 
sequencing quality as a comma separated string.

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

raw
:   specify how the sequencing quality should be returned. 
If `TRUE`, return the raw value. 
If `FALSE`, do WHAT?




Details
-------------------

Once a repertoire `data.frame` has been processed with [readFastqDb](readFastqDb.md) and 
contains the fields `quality_sequence` and `quality_sequence_num`,
`sequenceAlignmentQuality` can be used to retrieve the quality scores 
from the already present field `quality_sequence_num`, without requiring 
again the `fastq` file, and report them as a `data.frame` with sequencing 
qualities per position, not as a string. This is done setting `raw=TRUE`. 
This `data.frame` with qualities per position can be used to generate figures, 
for example.









