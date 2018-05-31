**padSeqEnds** - *Pads ragged ends of aligned DNA sequences*

Description
--------------------

`padSeqEnds` takes a vector of DNA sequences, as character strings,
and appends the ends of each sequence with an appropriate number of `"N"` 
characters to create a sequence vector with uniform lengths.


Usage
--------------------
```
padSeqEnds(seq, len = NULL, start = FALSE, pad_char = "N")
```

Arguments
-------------------

seq
:   character vector of DNA sequence strings.

len
:   length to pad to. Only applies if longer than the maximum length of
the data in `seq`.

start
:   if `TRUE` pad the beginning of each sequence instead of the end.

pad_char
:   character to use for padding.




Value
-------------------

A modified `seq` vector with padded sequences.



Examples
-------------------

```R
# Default behavior uniformly pads ragged ends
seq <- c("CCCCTGGG", "ACCCTG", "CCCC")
padSeqEnds(seq)

```


```
[1] "CCCCTGGG" "ACCCTGNN" "CCCCNNNN"

```


```R

# Pad to fixed length
padSeqEnds(seq, len=15)

```


```
[1] "CCCCTGGGNNNNNNN" "ACCCTGNNNNNNNNN" "CCCCNNNNNNNNNNN"

```


```R

# Add padding to the beginning of the sequences instead of the ends
padSeqEnds(seq, start=TRUE)

```


```
[1] "CCCCTGGG" "NNACCCTG" "NNNNCCCC"

```


```R
padSeqEnds(seq, len=15, start=TRUE)
```


```
[1] "NNNNNNNCCCCTGGG" "NNNNNNNNNACCCTG" "NNNNNNNNNNNCCCC"

```



See also
-------------------

See [maskSeqEnds](maskSeqEnds.md) for creating uniform masking from existing masking.



