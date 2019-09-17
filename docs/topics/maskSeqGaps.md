**maskSeqGaps** - *Masks gap characters in DNA sequences*

Description
--------------------

`maskSeqGaps` substitutes gap characters, `c("-", ".")`, with `"N"` 
in a vector of DNA sequences.


Usage
--------------------
```
maskSeqGaps(seq, mask_char = "N", outer_only = FALSE)
```

Arguments
-------------------

seq
:   character vector of DNA sequence strings.

mask_char
:   character to use for masking.

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
# Mask with Ns
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


```R

# Mask with dashes
maskSeqGaps(c("ATG-C", "CC..C"), mask_char="-")
```


```
[1] "ATG-C" "CC--C"

```



See also
-------------------

See [maskSeqEnds](maskSeqEnds.md) for masking ragged edges.






