





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

```


```


```


```R
print(x)
```


```
    PARENT   CHILD COUNT   EXPECTED PVALUE
1      IgA     IgA    39  62.200000    0.9
2      IgA IgA,IgG     3   3.300000    0.4
3      IgA     IgG     2   3.900000    0.8
4  IgA,IgG     IgA    29   8.400000    0.0
5  IgA,IgG IgA,IgG     1   1.500000    0.5
6  IgA,IgG     IgG    24   2.666667    0.0
7  IgD,IgG     IgG     8   2.166667    0.0
8      IgG     IgA     1   4.900000    1.0
9      IgG IgD,IgG     1   1.000000    0.0
10     IgG     IgG   112 135.100000    1.0

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.



