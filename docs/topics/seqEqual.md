





**seqEqual** - *Test DNA sequences for equality.*

Description
--------------------

`seqEqual` checks if two DNA sequences are identical.


Usage
--------------------
```
seqEqual(seq1, seq2, ignore = character())
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
seqEqual("ATG-C", "AT--C")

```


```
[1] TRUE

```


```R
seqEqual("ATGGC", "ATGGN")

```


```
[1] TRUE

```


```R
seqEqual("AT--T", "ATGGC")

```


```
[1] FALSE

```


```R

# Ignore only Ns
seqEqual("ATG-C", "AT--C", ignore="N")

```


```
[1] FALSE

```


```R
seqEqual("ATGGC", "ATGGN", ignore="N")

```


```
[1] TRUE

```


```R
seqEqual("AT--T", "ATGGC", ignore="N")
```


```
[1] FALSE

```



See also
-------------------

Used by [collapseDuplicates](collapseDuplicates.md).


