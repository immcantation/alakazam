





**maskSeqGaps** - *Masks gap characters in DNA sequences*

Description
--------------------

`maskSeqGaps` substitutes gap characters, `c("-", ".")`, with `"N"` 
in a vector of DNA sequences.


Usage
--------------------
```
maskSeqGaps(seq, outer_only = FALSE)
```

Arguments
-------------------

seq
:   a character vector of DNA sequence strings.

outer_only
:   if `TRUE` replace only contiguous leading and trailing gaps;
if `FALSE` replace all gap characters.




Value
-------------------

A modified `seq` vector with `"N"` in place of `c("-", ".")` 
characters.



Examples
-------------------

```R
maskSeqGaps(c("ATG-C", "CC..C"))

```


```
[1] "ATGNC" "CCNNC"

```


```R
maskSeqGaps("--ATG-C-")

```


```
[1] "NNATGNCN"

```


```R
maskSeqGaps("--ATG-C-", outer_only=TRUE)
```


```
[1] "NNATG-CN"

```



See also
-------------------

See [maskSeqEnds](maskSeqEnds.md) for masking ragged edges.



