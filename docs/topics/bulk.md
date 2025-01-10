**bulk** - *Calculates the average bulkiness of amino acid sequences*

Description
--------------------

`bulk` calculates the average bulkiness score of amino acid sequences. 
Non-informative positions are excluded, where non-informative is defined as any 
character in `c("X", "-", ".", "*")`.


Usage
--------------------
```
bulk(seq, bulkiness = NULL)
```

Arguments
-------------------

seq
:   vector of strings containing amino acid sequences.

bulkiness
:   named numerical vector defining bulkiness scores for 
each amino acid, where names are single-letter amino acid 
character codes. If `NULL`, then the Zimmerman et al, 1968
scale is used.




Value
-------------------

A vector of bulkiness scores for the sequence(s).


References
-------------------


1. Zimmerman JM, Eliezer N, Simha R. The characterization of amino acid sequences 
in proteins by statistical methods. J Theor Biol 21, 170-201 (1968).




Examples
-------------------

```R
# Default bulkiness scale
seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
bulk(seq)

```


```
[1] 14.46227 16.58857

```


```R

# Use the Grantham, 1974 side chain volumn scores from the seqinr package
library(seqinr)
data(aaindex)
x <- aaindex[["GRAR740103"]]$I
# Rename the score vector to use single-letter codes
names(x) <- translateStrings(names(x), ABBREV_AA)
# Calculate average volume
bulk(seq, bulkiness=x)

```


```
[1] 77.34091 93.71429

```



See also
-------------------

For additional size related indices see [aaindex](http://www.rdocumentation.org/packages/seqinr/topics/aaindex).






