**gravy** - *Calculates the hydrophobicity of amino acid sequences*

Description
--------------------

`gravy` calculates the Grand Average of Hydrophobicity (gravy) index 
of amino acid sequences using the method of Kyte & Doolittle. Non-informative
positions are excluded, where non-informative is defined as any character in 
`c("X", "-", ".", "*")`.


Usage
--------------------
```
gravy(seq, hydropathy = NULL)
```

Arguments
-------------------

seq
:   vector of strings containing amino acid sequences.

hydropathy
:   named numerical vector defining hydropathy index values for 
each amino acid, where names are single-letter amino acid 
character codes. If `NULL`, then the Kyte & Doolittle
scale is used.




Value
-------------------

A vector of gravy scores for the sequence(s).


References
-------------------


1. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character 
of a protein. J Mol Biol. 157, 105-32 (1982).




Examples
-------------------

```R
# Default scale
seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
gravy(seq)

```


```
[1] -0.9181818 -0.6571429

```


```R

# Use the Kidera et al, 1985 scores from the seqinr package
library(seqinr)
data(aaindex)
x <- aaindex[["KIDA850101"]]$I
# Rename the score vector to use single-letter codes
names(x) <- translateStrings(names(x), ABBREV_AA)
# Calculate hydrophobicity
gravy(seq, hydropathy=x)
```


```
[1] 0.3240909 0.3628571

```



See also
-------------------

For additional hydrophobicity indices see `[aaindex](http://www.rdocumentation.org/packages/seqinr/topics/aaindex)`.






