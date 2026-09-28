**seqMismatchPositions** - *Locate mismatches between sample and germline sequences*

Description
--------------------

`seqMismatchPositions` identifies Hamming-style mismatch positions between
paired sample and germline sequences, excluding ignored characters.


Usage
--------------------
```
seqMismatchPositions(samples, germlines, ignore = c("N", ".", "-"))
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
Default is to ignore `c("N", ".", "-")`.




Value
-------------------

List of integer vectors containing 1-based mismatch positions.


Details
-------------------

Comparisons are case-insensitive. Sequences of unequal length are
compared through the length of the shorter sequence.



Examples
-------------------

```R
seqMismatchPositions(c("ATGGC", "ATGGN", "TTGGG"), "ATGGC")

```


```
[[1]]
integer(0)

[[2]]
integer(0)

[[3]]
[1] 1 5


```



See also
-------------------

[seqMismatchCount](seqMismatchCount.md), [seqMismatchMatrix](seqMismatchMatrix.md)






