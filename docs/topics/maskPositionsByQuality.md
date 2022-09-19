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
fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")
db <- readFastqDb(db, fastq_file, quality_offset=-33)
maskPositionsByQuality(db, min_quality=90, quality_num="quality_alignment_num")
```

*Number of masked sequences: 1*
```
# A tibble: 1 × 37
  sequenc…¹ seque…² rev_c…³ produ…⁴ v_call d_call j_call seque…⁵ germl…⁶ junct…⁷
  <chr>     <chr>   <lgl>   <lgl>   <chr>  <chr>  <chr>  <chr>   <chr>   <chr>  
1 CGCTTTTC… GGCTTT… FALSE   TRUE    IGHV4… IGHD3… IGHJ4… CAGCTG… CAGCTG… TGTGCG…
# … with 27 more variables: junction_aa <chr>, v_cigar <chr>, d_cigar <chr>,
#   j_cigar <chr>, stop_codon <lgl>, vj_in_frame <lgl>, locus <chr>,
#   junction_length <int>, np1_length <int>, np2_length <int>,
#   v_sequence_start <int>, v_sequence_end <int>, v_germline_start <int>,
#   v_germline_end <int>, d_sequence_start <int>, d_sequence_end <int>,
#   d_germline_start <int>, d_germline_end <int>, j_sequence_start <int>,
#   j_sequence_end <int>, j_germline_start <int>, j_germline_end <int>, …

```



See also
-------------------

[readFastqDb](readFastqDb.md) and [getPositionQuality](getPositionQuality.md)






