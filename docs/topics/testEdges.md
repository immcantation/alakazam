**testEdges** - *Tests for parent-child annotation enchrichment in lineage trees*

Description
--------------------

`testEdges` performs a permutation test on a set of lineage trees to determine
the significance of an annotation's association with parent-child relationships.


Usage
--------------------
```
testEdges(
graphs,
field,
indirect = FALSE,
exclude = c("Germline", NA),
nperm = 200,
progress = FALSE
)
```

Arguments
-------------------

graphs
:   list of igraph objects with vertex annotations.

field
:   string defining the annotation field to permute.

indirect
:   if `FALSE` count direct connections (edges) only. If 
`TRUE` walk through any nodes with annotations specified in 
the `argument` to count indirect connections. Specifying
`indirect=TRUE` with `exclude=NULL` will have no effect.

exclude
:   vector of strings defining `field` values to exclude from 
permutation.

nperm
:   number of permutations to perform.

progress
:   if `TRUE` show a progress bar.




Value
-------------------

An [EdgeTest](EdgeTest-class.md) object containing the test results and permutation
realizations.



Examples
-------------------

```R
# Define example tree set
graphs <- ExampleTrees[1-10]

# Perform edge test on isotypes
x <- testEdges(graphs, "c_call", nperm=10)
print(x)
```


```
      parent     child count   expected pvalue
1       IGHA      IGHA    39  62.500000    0.9
2       IGHA IGHA,IGHG     3   3.500000    0.6
3       IGHA      IGHG     2   4.400000    0.7
4  IGHA,IGHG      IGHA    29   5.571429    0.0
5  IGHA,IGHG IGHA,IGHG     1   2.000000    1.0
6  IGHA,IGHG      IGHG    24   2.000000    0.0
7  IGHD,IGHG      IGHG     8   1.500000    0.0
8       IGHG      IGHA     1   4.900000    1.0
9       IGHG IGHD,IGHG     1   1.000000    0.0
10      IGHG      IGHG   112 135.500000    1.0

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.






