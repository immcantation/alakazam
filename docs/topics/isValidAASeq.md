





**isValidAASeq** - *Validate amino acid sequences*

Description
--------------------

`isValidAASeq` checks that a set of sequences are valid non-ambiguous 
amino acid sequences. A sequence is considered valid if it contains only 
characters in the the non-ambiguous IUPAC character set or any characters in 
`c("X", ".", "-", "*")`.

Usage
--------------------

```
isValidAASeq(seq)
```

Arguments
-------------------

seq
:   character vector of sequences to check.



Value
-------------------

A logical vector with `TRUE` for each valid amino acid sequences 
and `FALSE` for each invalid sequence.



Examples
-------------------

```R
seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVR--XX", "CARJ", "10") 
isValidAASeq(seq)
```


```
[1]  TRUE  TRUE FALSE FALSE

```



See also
-------------------

See [ABBREV_AA](ABBREV_AA.md) for the set of non-ambiguous amino acid characters.
See [IUPAC_AA](IUPAC_CODES.md) for the full set of ambiguous amino acid characters.



