





**aliphatic** - *Calculates the aliphatic index of amino acid sequences*

Description
--------------------

`aliphatic` calculates the aliphatic index of amino acid sequences using 
the method of Ikai. Non-informative positions are excluded, where non-informative 
is defined as any character in `c("X", "-", ".", "*")`.


Usage
--------------------
```
aliphatic(seq, normalize = TRUE)
```

Arguments
-------------------

seq
:   vector of strings containing amino acid sequences.

normalize
:   if `TRUE` then divide the aliphatic index of each amino acid 
sequence by the number of informative positions. Non-informative 
position are defined by the presence any character in 
`c("X", "-", ".", "*")`. If `FALSE` then return the raw
aliphatic index.



Value
-------------------

A vector of the aliphatic indices for the sequence(s).

References
-------------------


1. Ikai AJ. Thermostability and aliphatic index of globular proteins. 
J Biochem. 88, 1895-1898 (1980).




Examples
-------------------

```R
seq <- c("CARDRSTPWRRGIASTTVRTSW", NA, "XXTQMYVRT")
aliphatic(seq)
```


```
[1] 0.4000000        NA 0.4142857

```




