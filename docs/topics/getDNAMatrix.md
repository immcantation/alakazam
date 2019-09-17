**getDNAMatrix** - *Build a DNA distance matrix*

Description
--------------------

`getDNAMatrix` returns a Hamming distance matrix for IUPAC ambiguous
DNA characters with modifications for gap, `c("-", ".")`, and missing, 
`c("?")`, character values.


Usage
--------------------
```
getDNAMatrix(gap = -1)
```

Arguments
-------------------

gap
:   value to assign to characters in the set `c("-", ".")`.




Value
-------------------

A `matrix` of DNA character distances with row and column names 
indicating the character pair. By default, distances will be either 0 
(equivalent), 1 (non-equivalent or missing), or -1 (gap).



Examples
-------------------

```R
# Set gap characters to Inf distance
# Distinguishes gaps from Ns
getDNAMatrix()

```


```
   A  C  G  T  M  R  W  S  Y  K  V  H  D  B  N  -  . ?
A  0  1  1  1  0  0  0  1  1  1  0  0  0  1  0 -1 -1 1
C  1  0  1  1  0  1  1  0  0  1  0  0  1  0  0 -1 -1 1
G  1  1  0  1  1  0  1  0  1  0  0  1  0  0  0 -1 -1 1
T  1  1  1  0  1  1  0  1  0  0  1  0  0  0  0 -1 -1 1
M  0  0  1  1  0  0  0  0  0  1  0  0  0  0  0 -1 -1 1
R  0  1  0  1  0  0  0  0  1  0  0  0  0  0  0 -1 -1 1
W  0  1  1  0  0  0  0  1  0  0  0  0  0  0  0 -1 -1 1
S  1  0  0  1  0  0  1  0  0  0  0  0  0  0  0 -1 -1 1
Y  1  0  1  0  0  1  0  0  0  0  0  0  0  0  0 -1 -1 1
K  1  1  0  0  1  0  0  0  0  0  0  0  0  0  0 -1 -1 1
V  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 -1 -1 1
H  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0 -1 -1 1
D  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 -1 -1 1
B  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1 -1 1
N  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1 -1 1
- -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 1
. -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 1
?  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 0

```


```R

# Set gap characters to 0 distance
# Makes gap characters equivalent to Ns
getDNAMatrix(gap=0)
```


```
  A C G T M R W S Y K V H D B N - . ?
A 0 1 1 1 0 0 0 1 1 1 0 0 0 1 0 0 0 1
C 1 0 1 1 0 1 1 0 0 1 0 0 1 0 0 0 0 1
G 1 1 0 1 1 0 1 0 1 0 0 1 0 0 0 0 0 1
T 1 1 1 0 1 1 0 1 0 0 1 0 0 0 0 0 0 1
M 0 0 1 1 0 0 0 0 0 1 0 0 0 0 0 0 0 1
R 0 1 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1
W 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1
S 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 1
Y 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1
K 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1
V 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1
H 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
D 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
B 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
- 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
. 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
? 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0

```



See also
-------------------

Creates DNA distance matrix for [seqDist](seqDist.md).
See [getAAMatrix](getAAMatrix.md) for amino acid distances.






