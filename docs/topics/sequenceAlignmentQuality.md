**sequenceAlignmentQuality** - *Retrieve sequencing quality scores from a db file with `sequence_quality` information*

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
:   An AIRR data.frame

sequence_id
:   Name of the column in `data` that contains sequence 
identifiers to be matched to sequence identifiers in 
`fastq_file`.

sequence
:   Name of the column in `data` that contains sequence data.

sequence_alignment
:   Name of the column in `data` that contains 
IMGT aligned sequence data.

v_cigar
:   Name of the column in `data` that contains CIGAR 
strings for the V gene alignments.

d_cigar
:   Name of the column in `data` that contains CIGAR 
strings for the D gene alignments.

j_cigar
:   Name of the column in `data` that contains CIGAR 
strings for the J gene alignments.

np1_length
:   Name of the column in `data` that contains the number
of nucleotides between the V gene and first D gene 
alignments or between the V gene and J gene alignments.

np2_length
:   Name of the column in `data` that contains the number
of nucleotides between either the first D gene and J 
gene alignments or the first D gene and second D gene
alignments.

v_sequence_end
:   Name of the column in `data` that contains the 
end position of the V gene in `sequence`.

d_sequence_end
:   Name of the column in `data` that contains the 
end position of the D gene in `sequence`.

raw
:   Select how the sequencing quality should be returned. 
`num` to store the quality scores as strings of 
comma separated numeric values. Use `phred` to have
the function return the scores as phred (ascii) scores. 
Use `both` to retrieve both.

sequence_quality
:   Name of the column in `data` that contains 
sequencing quality (phred scores).

sequence_quality_num
:   Name of the column in `data` that contains 
sequencing quality, as comma separated string.




Details
-------------------

Once a repertoire `data.frame` has been processed with `readFasqDb` and 
contains the fields `quality_sequence` and `quality_sequence_num`,
`sequenceAlignmentQuality` can
be used to retrieve the quality scores from the already present field
 `sequence_quality`, without requiring again the `fastq` file, 
and report them as a `data.frame` with sequencing qualities per position, 
not as a string. This is done setting `raw=TRUE`. This `data.frame` 
with qualities per position can be used to generate figures, for example.









