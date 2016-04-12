





**maskSeqEnds** - *Masks ragged leading and trailing edges of aligned DNA sequences*

Description
--------------------

`maskSeqEnds` takes a vector of DNA sequences, as character strings,
and replaces the leading and trailing characters with `"N"` characters to create 
a sequence vector with uniformly masked outer sequence segments.


Usage
--------------------
```
maskSeqEnds(seq, max_mask = NULL, trim = FALSE)
```

Arguments
-------------------

seq
:   a character vector of DNA sequence strings.

max_mask
:   the maximum number of characters to mask. If set to 0 then
no masking will be performed. If set to `NULL` then the upper 
masking bound will be automatically determined from the maximum 
number of observed leading or trailing `"N"` characters amongst 
all strings in `seq`.

trim
:   if `TRUE` leading and trailing characters will be cut rather 
than masked with `"N"` characters.



Value
-------------------

A modified `seq` vector with masked (or optionally trimmed) sequences.



Examples
-------------------

```R
# Default behavior uniformly masks ragged ends
seq <- c("CCCCTGGG", "NAACTGGN", "NNNCTGNN")
maskSeqEnds(seq)

```


```
[1] "NNNCTGNN" "NNNCTGNN" "NNNCTGNN"

```


```R

# Does nothing
maskSeqEnds(seq, max_mask=0)

```


```
[1] "CCCCTGGG" "NAACTGGN" "NNNCTGNN"

```


```R

# Cut ragged sequence ends
maskSeqEnds(seq, trim=TRUE)

```


```
[1] "CTG" "CTG" "CTG"

```


```R

# Set max_mask to limit extent of masking and trimming
maskSeqEnds(seq, max_mask=1)

```


```
[1] "NCCCTGGN" "NAACTGGN" "NNNCTGNN"

```


```R
maskSeqEnds(seq, max_mask=1, trim=TRUE)
```


```
[1] "CCCTGG" "AACTGG" "NNCTGN"

```



See also
-------------------

See [maskSeqGaps](maskSeqGaps.md) for masking internal gaps.



