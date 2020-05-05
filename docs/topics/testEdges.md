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
      parent     child count   expected    pvalue
1       IGHA      IGHA    39  58.000000 0.9000000
2       IGHA IGHA,IGHG     3   3.222222 0.3333333
3       IGHA      IGHG     2   6.700000 0.8000000
4  IGHA,IGHG      IGHA    29  10.000000 0.0000000
5  IGHA,IGHG IGHA,IGHG     1   2.000000 1.0000000
6  IGHA,IGHG      IGHG    24   2.250000 0.0000000
7  IGHD,IGHG      IGHG     8   2.750000 0.0000000
8       IGHG      IGHA     1   6.000000 0.8000000
9       IGHG IGHD,IGHG     1   1.000000 0.0000000
10      IGHG      IGHG   112 134.900000 1.0000000

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.






