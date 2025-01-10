**polar** - *Calculates the average polarity of amino acid sequences*

Description
--------------------

`polar` calculates the average polarity score of amino acid sequences. 
Non-informative positions are excluded, where non-informative is defined as any 
character in `c("X", "-", ".", "*")`.


Usage
--------------------
```
polar(seq, polarity = NULL)
```

Arguments
-------------------

seq
:   vector of strings containing amino acid sequences.

polarity
:   named numerical vector defining polarity scores for 
each amino acid, where names are single-letter amino acid 
character codes. If `NULL`, then the Grantham, 1974
scale is used.




Value
-------------------

A vector of bulkiness scores for the sequence(s).


References
-------------------


1. Grantham R. Amino acid difference formula to help explain protein evolution. 
Science 185, 862-864 (1974).




Examples
-------------------

```R
# Default scale
seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
polar(seq)

```


```
[1] 8.55 8.00

```


```R

# Use the Zimmerman et al, 1968 polarity scale from the seqinr package
library(seqinr)
data(aaindex)
x <- aaindex[["ZIMJ680103"]]$I
# Rename the score vector to use single-letter codes
names(x) <- translateStrings(names(x), ABBREV_AA)
# Calculate polarity
polar(seq, polarity=x)

```


```
[1] 14.94864  8.86000

```



See also
-------------------

For additional size related indices see `[aaindex](http://www.rdocumentation.org/packages/seqinr/topics/aaindex)`.






