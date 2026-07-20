**fastDist** - *Faster calculation of pairwise distances between sequences of the same length and contain only "ACTGN?"*

Description
--------------------

`fastDist` calculates all pairwise distance between a set of sequences of the same length and contain only "ACTGN?".


Usage
--------------------
```
fastDist(seqs)
```

Arguments
-------------------

seqs
:   character vector containing a DNA sequences.




Value
-------------------

Packed lower triangular matrix of distance between each entry in `seq`. 
If `seq` is a named vector, row and columns names will be added 
accordingly.



Examples
-------------------

```R
fastDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT?NC", E="NTGG?"))

```


```
  A B C D
B 1      
C 1 0    
D 1 2 2  
E 1 1 1 2

```








