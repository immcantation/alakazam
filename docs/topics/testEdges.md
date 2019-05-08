**testEdges** - *Tests for parent-child annotation enchrichment in lineage trees*

Description
--------------------

`testEdges` performs a permutation test on a set of lineage trees to determine
the significance of an annotation's association with parent-child relationships.


Usage
--------------------
```
testEdges(graphs, field, indirect = FALSE, exclude = c("Germline", NA),
nperm = 200, progress = FALSE)
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
x <- testEdges(graphs, "ISOTYPE", nperm=10)
print(x)
```


```
    PARENT   CHILD COUNT   EXPECTED PVALUE
1      IgA     IgA    39  63.300000    1.0
2      IgA IgA,IgG     3   4.000000    0.7
3      IgA     IgG     2   6.100000    1.0
4  IgA,IgG     IgA    29   4.714286    0.0
5  IgA,IgG IgA,IgG     1   3.000000    1.0
6  IgA,IgG     IgG    24   1.800000    0.0
7  IgD,IgG     IgG     8   1.666667    0.0
8      IgG     IgA     1   3.600000    1.0
9      IgG IgD,IgG     1   1.000000    0.0
10     IgG     IgG   112 135.200000    1.0

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.



