





**countPatterns** - *Count sequence patterns*

Description
--------------------

`countPatterns` counts the fraction of times a set of character patterns occur 
in a set of sequences.

Usage
--------------------

```
countPatterns(seq, patterns, nt = FALSE, trim = FALSE, label = "REGION")
```

Arguments
-------------------

seq
:   character vector of either DNA or amino acid sequences.

patterns
:   list of sequence patterns to count in each sequence. If the 
list is named, then names will be assigned as the column names of 
output data.frame.

nt
:   if `TRUE` then `seq` are DNA sequences and and will be 
translated before performing the pattern search.

trim
:   if `TRUE` remove the first and last codon or amino acid from 
each sequence before the pattern search. If `FALSE` do
not modify the input sequences.

label
:   string defining a label to add as a prefix to the output 
column names.



Value
-------------------

A data.frame containing the fraction of times each sequence pattern was 
found.



Examples
-------------------

```R
seq <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
"TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
"TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
patterns <- c("A", "V", "[LI]")
names(patterns) <- c("ARG", "VAL", "ISO_LEU")
countPatterns(seq, patterns, nt=TRUE, trim=TRUE, label="CDR3")
```


```
   CDR3_ARG CDR3_VAL CDR3_ISO_LEU
1 0.1250000        0    0.0000000
2 0.1111111        0    0.1111111
3 0.1111111        0    0.0000000

```




