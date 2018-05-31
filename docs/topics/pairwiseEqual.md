**pairwiseEqual** - *Calculate pairwise equivalence between sequences*

Description
--------------------

`pairwiseEqual` determined pairwise equivalence between a pairs in a 
set of sequences, excluding ambiguous positions (Ns and gaps).


Usage
--------------------
```
pairwiseEqual(seq)
```

Arguments
-------------------

seq
:   character vector containing a DNA sequences.




Value
-------------------

A logical matrix of equivalence between each entry in `seq`. 
Values are `TRUE` when sequences are equivalent and `FALSE`
when they are not.



Examples
-------------------

```R
# Gaps and Ns will match any character
seq <- c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C", E="NTGGG")
d <- pairwiseEqual(seq)
rownames(d) <- colnames(d) <- seq
d
```


```
      ATGGC ATGGG ATGGG AT--C NTGGG
ATGGC  TRUE FALSE FALSE  TRUE FALSE
ATGGG FALSE  TRUE  TRUE FALSE  TRUE
ATGGG FALSE  TRUE  TRUE FALSE  TRUE
AT--C  TRUE FALSE FALSE  TRUE FALSE
NTGGG FALSE  TRUE  TRUE FALSE  TRUE

```



See also
-------------------

Uses [seqEqual](seqEqual.md) for testing equivalence between pairs.
See [pairwiseDist](pairwiseDist.md) for generating a sequence distance matrix.



