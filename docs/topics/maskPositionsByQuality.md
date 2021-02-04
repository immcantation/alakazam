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
fastq_file <- system.file("extdata", "test_seq.fastq", package="alakazam")
db <- readFastqDb(db, fastq_file, quality_offset=-33)
maskPositionsByQuality(db, min_quality=90, quality_num="quality_alignment_num")
```

*Number of masked sequences: 1*
```
[90m# A tibble: 1 x 37[39m
  sequence_id sequence rev_comp productive v_call d_call j_call sequence_alignm…
  [3m[90m<chr>[39m[23m       [3m[90m<chr>[39m[23m    [3m[90m<lgl>[39m[23m    [3m[90m<lgl>[39m[23m      [3m[90m<chr>[39m[23m  [3m[90m<chr>[39m[23m  [3m[90m<chr>[39m[23m  [3m[90m<chr>[39m[23m           
[90m1[39m CGCTTTTCGG… GGCTTTC… FALSE    TRUE       IGHV4… IGHD3… IGHJ4… CAGCTGCACCTGCAG…
[90m# … with 29 more variables: germline_alignment [3m[90m<chr>[90m[23m, junction [3m[90m<chr>[90m[23m,[39m
[90m#   junction_aa [3m[90m<chr>[90m[23m, v_cigar [3m[90m<chr>[90m[23m, d_cigar [3m[90m<chr>[90m[23m, j_cigar [3m[90m<chr>[90m[23m,[39m
[90m#   stop_codon [3m[90m<lgl>[90m[23m, vj_in_frame [3m[90m<lgl>[90m[23m, locus [3m[90m<chr>[90m[23m, junction_length [3m[90m<int>[90m[23m,[39m
[90m#   np1_length [3m[90m<int>[90m[23m, np2_length [3m[90m<int>[90m[23m, v_sequence_start [3m[90m<int>[90m[23m,[39m
[90m#   v_sequence_end [3m[90m<int>[90m[23m, v_germline_start [3m[90m<int>[90m[23m, v_germline_end [3m[90m<int>[90m[23m,[39m
[90m#   d_sequence_start [3m[90m<int>[90m[23m, d_sequence_end [3m[90m<int>[90m[23m, d_germline_start [3m[90m<int>[90m[23m,[39m
[90m#   d_germline_end [3m[90m<int>[90m[23m, j_sequence_start [3m[90m<int>[90m[23m, j_sequence_end [3m[90m<int>[90m[23m,[39m
[90m#   j_germline_start [3m[90m<int>[90m[23m, j_germline_end [3m[90m<int>[90m[23m, consensus_count [3m[90m<int>[90m[23m,[39m
[90m#   duplicate_count [3m[90m<int>[90m[23m, prcons [3m[90m<chr>[90m[23m, quality_alignment_num [3m[90m<chr>[90m[23m,[39m
[90m#   sequence_alignment_masked [3m[90m<chr>[90m[23m[39m

```



See also
-------------------

[readFastqDb](readFastqDb.md) and [getPositionQuality](getPositionQuality.md)






