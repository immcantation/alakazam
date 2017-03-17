





**charge** - *Calculates the net charge of amino acid sequences.*

Description
--------------------

`charge` calculates the net charge of amino acid sequences using 
the method of Moore, 1985, with exclusion of the C-terminus and N-terminus charges.


Usage
--------------------
```
charge(seq, pH = 7.4, pK = NULL, normalize = FALSE)
```

Arguments
-------------------

seq
:   vector strings defining of amino acid sequences.

pH
:   environmental pH.

pK
:   named vector defining pK values for each charged amino acid,
where names are the single-letter amino acid character codes
`c("R", "H", "K", "D", "E", "C", "Y")`). If `NULL`, 
then the EMBOSS scale is used.

normalize
:   if `TRUE` then divide the net charge of each amino acid 
sequence by the number of informative positions. Non-informative 
position are defined by the presence any character in 
`c("X", "-", ".", "*")`. If `FALSE` then return the raw
net charge.




Value
-------------------

A vector of net charges for the sequence(s).


References
-------------------


1. Moore DS. Amino acid and peptide net charges: A simple calculational procedure. 
Biochem Educ. 13, 10-11 (1985).
1. [http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html](http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html)




Examples
-------------------

```R
seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT") 
# Unnormalized charge
charge(seq)

```


```
[1] 3.9266889 0.9980008

```


```R
# Normalized charge
charge(seq, normalize=TRUE)

```


```
[1] 0.1784859 0.1425715

```


```R

# Use the Murray et al, 2006 scores from the seqinr package
library(seqinr)
data(pK)
x <- setNames(pK[["Murray"]], rownames(pK))
# Calculate charge
charge(seq, pK=x)
```


```
[1] 3.8946562 0.9977872

```



See also
-------------------

For additional pK scales see `[pK](http://www.inside-r.org/packages/cran/seqinr/docs/pK)`.



