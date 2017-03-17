





**seqDist** - *Calculate distance between two sequences*

Description
--------------------

`seqDist` calculates the distance between two DNA sequences.


Usage
--------------------
```
seqDist(seq1, seq2, dist_mat = getDNAMatrix())
```

Arguments
-------------------

seq1
:   character string containing a DNA sequence.

seq2
:   character string containing a DNA sequence.

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

Numerical distance between `seq1` and `seq2`.



Examples
-------------------

```R
# Ungapped examples
seqDist("ATGGC", "ATGGG")

```


```
[1] 1

```


```R
seqDist("ATGGC", "ATG??")

```


```
[1] 2

```


```R

# Gaps will be treated as Ns with a gap=0 distance matrix
seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=0))

```


```
[1] 0

```


```R

# Gaps will be treated as universally non-matching characters with gap=1
seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=1))

```


```
[1] 2

```


```R

# Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))

```


```
[1] 1

```


```R

# Gaps of equivalent run lengths are not counted as gaps
seqDist("ATG-C", "ATG-C", dist_mat=getDNAMatrix(gap=-1))

```


```
[1] 0

```


```R

# Overlapping runs of gap characters are counted as a single gap
seqDist("ATG-C", "AT--C", dist_mat=getDNAMatrix(gap=-1))

```


```
[1] 1

```


```R
seqDist("A-GGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))

```


```
[1] 1

```


```R
seqDist("AT--C", "AT--C", dist_mat=getDNAMatrix(gap=-1))

```


```
[1] 0

```


```R

# Discontiguous runs of gap characters each count as separate gaps
seqDist("-TGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))
```


```
[1] 2

```



See also
-------------------

Nucleotide distance matrix may be built with 
[getDNAMatrix](getDNAMatrix.md). Amino acid distance matrix may be built
with [getAAMatrix](getAAMatrix.md). Used by [pairwiseDist](pairwiseDist.md) for generating
distance matrices. See [seqEqual](seqEqual.md) for testing sequence equivalence.



