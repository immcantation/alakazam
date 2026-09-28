**seqMismatchCount** - *Count or locate mismatches between sample and germline sequences*

Description
--------------------

`seqMismatchCount` counts Hamming-style mismatches between paired sample
and germline sequences, `seqMismatchMatrix` counts them between every
sample and every germline, and `seqMismatchPositions` returns the
mismatch positions of paired sequences.


Usage
--------------------
```
seqMismatchCount(
samples,
germlines,
ignore = c("N", "-", ".", "?"),
count_trailing = FALSE
)
```
```
seqMismatchMatrix(
samples,
germlines,
ignore = c("N", "-", ".", "?"),
count_trailing = FALSE
)
```
```
seqMismatchPositions(
samples,
germlines,
ignore = c("N", "-", ".", "?"),
count_trailing = FALSE
)
```

Arguments
-------------------

samples
:   character vector of sample sequences.

germlines
:   character vector of germline sequences. For
`seqMismatchCount` and `seqMismatchPositions`,
a single germline is recycled across all samples.

ignore
:   vector of characters to ignore, in either sequence.
Default is to ignore `c("N", "-", ".", "?")`.

count_trailing
:   if `TRUE`, sample positions past the end of a
shorter germline count as mismatches, so a germline
gains nothing from ending early. If `FALSE`,
sequences are compared only through the length of the
shorter one.




Value
-------------------

`seqMismatchCount`: an integer vector of mismatch counts.
`seqMismatchMatrix`: an integer matrix of mismatch counts, samples
in rows and germlines in columns.
`seqMismatchPositions`: a list of integer vectors of 1-based
mismatch positions.


Details
-------------------

Comparisons are case-insensitive. A missing (`NA`) sample or
germline gives `NA`.



Examples
-------------------

```R
seqMismatchCount(c("ATGGC", "ATGGN"), "ATGGC")

```


```
[1] 0 0

```


```R
seqMismatchMatrix(c("ATGGC", "ATGGN"), c("ATGGC", "ATGGG"))

```


```
     [,1] [,2]
[1,]    0    1
[2,]    0    0

```


```R
seqMismatchPositions("ATGGCA", "ATGG")

```


```
[[1]]
integer(0)


```


```R

# A germline that ends early is not rewarded for it
seqMismatchMatrix("ATGGCA", c(full="ATGGCC", short="ATGG"))

```


```
     full short
[1,]    1     0

```


```R
seqMismatchMatrix("ATGGCA", c(full="ATGGCC", short="ATGG"), count_trailing=TRUE)

```


```
     full short
[1,]    1     2

```








