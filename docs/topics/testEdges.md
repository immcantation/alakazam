





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
1      IgA     IgA    39  62.000000  1.000
2      IgA IgA,IgG     3   3.600000  0.600
3      IgA     IgG     2   4.375000  0.875
4  IgA,IgG     IgA    29   5.857143  0.000
5  IgA,IgG IgA,IgG     1   2.000000  1.000
6  IgA,IgG     IgG    24   4.750000  0.000
7  IgD,IgG     IgG     8   2.400000  0.000
8      IgG     IgA     1   4.700000  0.900
9      IgG IgD,IgG     1   1.000000  0.000
10     IgG     IgG   112 134.000000  1.000

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.



