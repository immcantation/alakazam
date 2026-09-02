**fastDistAA** - *Faster calculation of pairwise distances between amino acid sequences of the same length*

Description
--------------------

`fastDistAA` calculates all pairwise distances among a set of amino acid sequences of the same length. 
Amino acid sequences may contain the 20 standard amino acid characters and four special characters: (X, ., -, *).  
The characters `X`, `-`, and `.` match any character, whereas standard amino acids and stop codon `*` match only themselves.


Usage
--------------------
```
fastDistAA(seqs)
```

Arguments
-------------------

seqs
:   character vector containing an amino acid sequences.




Value
-------------------

Packed lower triangular matrix of distance between each entry in `seq`. 
If `seq` is a named vector, row and columns names will be added 
accordingly.



Examples
-------------------

```R
fastDistAA(c(A="AEHGC*X", B="AEHGGIC", C="AXGGGIC", D="ATTNC-E", E="N.TGG**"))

```


```
  A B C D
B 2      
C 3 1    
D 3 5 4  
E 3 4 4 4

```








