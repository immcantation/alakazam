**seqMismatchMatrixRcpp** - *Count mismatches between samples and germlines.*

Description
--------------------

`seqMismatchMatrixRcpp` counts Hamming-style mismatches between each sample
and each germline sequence, excluding ignored characters.


Usage
--------------------
```
seqMismatchMatrixRcpp(
samples,
germlines,
ignore = as.character(c("N", ".", "-"))
)
```

Arguments
-------------------

samples
:   character vector containing sample sequences.

germlines
:   character vector containing germline sequences.

ignore
:   vector of characters to ignore when counting mismatches.
Default is to ignore c("N", ".", "-").




Value
-------------------

Integer matrix of mismatch counts, with rows corresponding to
samples and columns corresponding to germlines.


Details
-------------------

Comparisons are case-insensitive. Sequences of unequal length are
compared through the length of the shorter sequence.









