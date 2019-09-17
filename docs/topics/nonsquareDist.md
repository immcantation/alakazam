**nonsquareDist** - *Calculate pairwise distances between sequences*

Description
--------------------

`nonsquareDist` calculates all pairwise distance between a set of sequences and a subset of it.


Usage
--------------------
```
nonsquareDist(seq, indx, dist_mat = getDNAMatrix())
```

Arguments
-------------------

seq
:   character vector containing a DNA sequences. The sequence vector needs to
be named.

indx
:   numeric vector contating the indices (a subset of indices of `seq`).

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

A matrix of numerical distance between each entry in `seq` and 
sequences specified by `indx` indices. 

Note that the input subsampled indices will be ordered ascendingly. Therefore, 
it is necassary to assign unique names to the input sequences, `seq`, 
to recover the input order later. Row and columns names will be added accordingly.

Amino acid distance matrix may be built with [getAAMatrix](getAAMatrix.md). 
Uses [seqDist](seqDist.md) for calculating distances between pairs.
See [pairwiseEqual](pairwiseEqual.md) for generating an equivalence matrix.



Examples
-------------------

```R
# Gaps will be treated as Ns with a gap=0 distance matrix
seq <- c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C")
pairwiseDist(seq, 
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

nonsquareDist(seq, indx=c(1,3), 
dist_mat=getDNAMatrix(gap=0))
```


```
  A B C D
A 0 1 1 0
C 1 0 0 1

```








