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
db <- airr::read_rearrangement(system.file("extdata", "example_quality.tsv", package="alakazam"))

```


```
[1mindexed[0m [32m0B[0m in [36m 0s[0m, [32m0B/s[0m[1mindexed[0m [32m1.00TB[0m in [36m 0s[0m, [32m960.67TB/s[0m                                                                                                                                                                                                                                [1mindexed[0m [32m0B[0m in [36m 0s[0m, [32m0B/s[0m[1mindexed[0m [32m1.00TB[0m in [36m 0s[0m, [32m1.02PB/s[0m                                                                                                                                                                                                                                
```


```R
fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")
db <- readFastqDb(db, fastq_file, quality_offset=-33)
maskPositionsByQuality(db, min_quality=90, quality_num="quality_alignment_num")
```

*Number of masked sequences: 1*
```
[38;5;246m# A tibble: 1 × 37[39m
  sequence_id sequence rev_comp productive v_call d_call j_call sequence_alignm… germline_alignm… junction junction_aa v_cigar d_cigar j_cigar stop_codon vj_in_frame locus junction_length np1_length np2_length v_sequence_start
  [3m[38;5;246m<chr>[39m[23m       [3m[38;5;246m<chr>[39m[23m    [3m[38;5;246m<lgl>[39m[23m    [3m[38;5;246m<lgl>[39m[23m      [3m[38;5;246m<chr>[39m[23m  [3m[38;5;246m<chr>[39m[23m  [3m[38;5;246m<chr>[39m[23m  [3m[38;5;246m<chr>[39m[23m            [3m[38;5;246m<chr>[39m[23m            [3m[38;5;246m<chr>[39m[23m    [3m[38;5;246m<chr>[39m[23m       [3m[38;5;246m<chr>[39m[23m   [3m[38;5;246m<chr>[39m[23m   [3m[38;5;246m<chr>[39m[23m   [3m[38;5;246m<lgl>[39m[23m      [3m[38;5;246m<lgl>[39m[23m       [3m[38;5;246m<chr>[39m[23m           [3m[38;5;246m<int>[39m[23m      [3m[38;5;246m<int>[39m[23m      [3m[38;5;246m<int>[39m[23m            [3m[38;5;246m<int>[39m[23m
[38;5;250m1[39m CGCTTTTCGG… GGCTTTC… FALSE    TRUE       IGHV4… IGHD3… IGHJ4… CAGCTGCACCTGCAG… CAGCTGCAGCTGCAG… TGTGCGA… CARGTDLVTG… 93S8=1… 403S9N… 425S4N… FALSE      TRUE        IGH                54         10          5               94
[38;5;246m# … with 16 more variables: v_sequence_end <int>, v_germline_start <int>, v_germline_end <int>, d_sequence_start <int>, d_sequence_end <int>, d_germline_start <int>, d_germline_end <int>, j_sequence_start <int>,[39m
[38;5;246m#   j_sequence_end <int>, j_germline_start <int>, j_germline_end <int>, consensus_count <int>, duplicate_count <int>, c_call <chr>, quality_alignment_num <chr>, sequence_alignment_masked <chr>[39m

```



See also
-------------------

[readFastqDb](readFastqDb.md) and [getPositionQuality](getPositionQuality.md)






