





**pairwiseDist** - *Calculate pairwise distances between sequences*

Description
--------------------

`pairwiseDist` calculates all pairwise distance between a set of sequences.


Usage
--------------------
```
pairwiseDist(seq, dist_mat = getDNAMatrix())
```

Arguments
-------------------

seq
:   character vector containing a DNA sequences.

dist_mat
:   Character distance matrix. Defaults to a Hamming distance 
matrix returned by [getDNAMatrix](getDNAMatrix.md). If gap 
characters, `c("-", ".")`, are assigned a value of -1 
in `dist_mat` then contiguous gaps of any run length,
which are not present in both sequences, will be counted as a 
distance of 1. Meaning, indels of any length will increase
the sequence distance by 1. Gap values other than -1 will 
return a distance that does not consider indels as a special case.



Value
-------------------

A matrix of numerical distance between each entry in `seq`. 
If `seq` is a named vector, row and columns names will be added 
accordingly.



Examples
-------------------

```R
# Gaps will be treated as Ns with a gap=0 distance matrix
pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
dist_mat=getDNAMatrix(gap=0))

```


```
  A B C D
A 0 1 1 0
B 1 0 0 1
C 1 0 0 1
D 0 1 1 0

```


```R

# Gaps will be treated as universally non-matching characters with gap=1
pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
dist_mat=getDNAMatrix(gap=1))

```


```
  A B C D
A 0 1 1 2
B 1 0 0 3
C 1 0 0 3
D 2 3 3 0

```


```R

# Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
dist_mat=getDNAMatrix(gap=-1))
```


```
  A B C D
A 0 1 1 1
B 1 0 0 2
C 1 0 0 2
D 1 2 2 0

```



See also
-------------------

Nucleotide distance matrix may be built with [getDNAMatrix](getDNAMatrix.md). 
Amino acid distance matrix may be built with [getAAMatrix](getAAMatrix.md). 
Uses [seqDist](seqDist.md) for calculating distances between pairs.
See [pairwiseEqual](pairwiseEqual.md) for generating an equivalence matrix.



