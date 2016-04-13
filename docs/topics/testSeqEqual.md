





**testSeqEqual** - *Test DNA sequences for equality.*

Description
--------------------

`testSeqEqual` checks if two DNA sequences are identical.


Usage
--------------------
```
testSeqEqual(seq1, seq2, ignore = c("N", "-", ".", "?"), rcpp = F)
```

Arguments
-------------------

seq1
:   character string containing a DNA sequence.

seq2
:   character string containing a DNA sequence.

ignore
:   vector of characters to ignore when testing for equality.

rcpp
:   use the Rcpp version of the code



Value
-------------------

Returns `TRUE` if sequences are equal and `FALSE` if they are not.
Sequences of unequal length will always return `FALSE` regardless of
their character values.



Examples
-------------------

```R
# Ignore gaps
testSeqEqual("ATG-C", "AT--C")

```


```
[1] TRUE

```


```R
testSeqEqual("ATGGC", "ATGGN")

```


```
[1] TRUE

```


```R
testSeqEqual("AT--T", "ATGGC")

```


```
[1] FALSE

```


```R

# Ignore only Ns
testSeqEqual("ATG-C", "AT--C", ignore="N")

```


```
[1] FALSE

```


```R
testSeqEqual("ATGGC", "ATGGN", ignore="N")

```


```
[1] TRUE

```


```R
testSeqEqual("AT--T", "ATGGC", ignore="N")
```


```
[1] FALSE

```



See also
-------------------

Used by [collapseDuplicates](collapseDuplicates.md).



