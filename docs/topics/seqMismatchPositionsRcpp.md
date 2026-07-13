**seqMismatchPositionsRcpp** - *Locate mismatches between sample and germline sequences.*

Description
--------------------

`seqMismatchPositionsRcpp` identifies Hamming-style mismatch positions between
paired sample and germline sequences, excluding ignored characters.


Usage
--------------------
```
seqMismatchPositionsRcpp(
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
:   character vector containing germline sequences. If length
one, the germline is recycled across all samples.

ignore
:   vector of characters to ignore when locating mismatches.
Default is to ignore c("N", ".", "-").




Value
-------------------

List of integer vectors containing 1-based mismatch positions.


Details
-------------------

Comparisons are case-insensitive. Sequences of unequal length are
compared through the length of the shorter sequence.









