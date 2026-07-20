**seqMismatchCountRcpp** - *Count mismatches between sample and germline sequences.*

Description
--------------------

`seqMismatchCountRcpp` counts Hamming-style mismatches between paired sample and
germline sequences, excluding ignored characters.


Usage
--------------------
```
seqMismatchCountRcpp(
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
:   vector of characters to ignore when counting mismatches.
Default is to ignore c("N", ".", "-").




Value
-------------------

Integer vector of mismatch counts.


Details
-------------------

Comparisons are case-insensitive. Sequences of unequal length are
compared through the length of the shorter sequence.



Examples
-------------------

```R
# Single germline recycled across samples
seqMismatchCountRcpp(c("ATGGC", "ATGGN"), "ATGGC")

```


```
[1] 0 0

```


```R

# Paired germlines and custom ignore characters
seqMismatchCountRcpp(c("ATG-C", "AT--C"), c("ATGGC", "ATGGG"), ignore="N")

```


```
[1] 1 3

```








