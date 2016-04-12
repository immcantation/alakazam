





**rcpp_testSeqEqual** - *Test DNA sequences for equality.*

Description
--------------------

`rcpp_testSeqEqual` checks if two DNA sequences are identical.


Usage
--------------------
```
rcpp_testSeqEqual(seq1, seq2, ignore = character())
```

Arguments
-------------------

seq1
:   character string containing a DNA sequence.

seq2
:   character string containing a DNA sequence.

ignore
:   vector of characters to ignore when testing for equality.



Value
-------------------

Returns `TRUE` if sequences are equal and `FALSE` if they are not.
Sequences of unequal length will always return `FALSE` regardless of
their character values.



Examples
-------------------

```R
# Ignore gaps
rcpp_testSeqEqual("ATG-C", "AT--C")

```


```
[1] TRUE

```


```R
rcpp_testSeqEqual("ATGGC", "ATGGN")

```


```
[1] TRUE

```


```R
rcpp_testSeqEqual("AT--T", "ATGGC")

```


```
[1] FALSE

```


```R

# Ignore only Ns
rcpp_testSeqEqual("ATG-C", "AT--C", ignore="N")

```


```
[1] FALSE

```


```R
rcpp_testSeqEqual("ATGGC", "ATGGN", ignore="N")

```


```
[1] TRUE

```


```R
rcpp_testSeqEqual("AT--T", "ATGGC", ignore="N")
```


```
[1] FALSE

```



See also
-------------------

Used by [collapseDuplicates](collapseDuplicates.md).



