**testEdges** - *Tests for parent-child annotation enrichment in lineage trees*

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
graphs <- ExampleTrees[1:10]

# Perform edge test on isotypes
x <- testEdges(graphs, "c_call", nperm=10)
print(x)

```


```
  parent child count expected pvalue
1   IGHA  IGHA     6        6      0
2   IGHG  IGHG     6        6      0

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.






