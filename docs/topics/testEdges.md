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
the `exclude` argument to count indirect connections. Specifying
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
1      IgA     IgA    39  57.900000    1.0
2      IgA IgA,IgG     3   2.900000    0.5
3      IgA     IgG     2   5.300000    0.9
4  IgA,IgG     IgA    29  10.000000    0.0
5  IgA,IgG IgA,IgG     1   2.000000    0.5
6  IgA,IgG     IgG    24   1.777778    0.0
7  IgD,IgG     IgG     8   4.500000    0.0
8      IgG     IgA     1   3.600000    0.9
9      IgG IgD,IgG     1   1.000000    0.0
10     IgG     IgG   112 134.100000    1.0

```


```R

# Perform edge test with inclusion of indirect connections
# Walks through germline, inferred, and nodes annotated as "IgM"
x <- testEdges(graphs, "ISOTYPE", nperm=10, indirect=TRUE, 
exclude=c("Germline", "IgM", NA))
print(x)
```


```
    PARENT   CHILD COUNT   EXPECTED PVALUE
1      IgA     IgA    41  70.900000    1.0
2      IgA IgA,IgG     4   4.200000    0.3
3      IgA     IgG     4   3.900000    0.3
4  IgA,IgG     IgA    31   1.857143    0.0
5  IgA,IgG IgA,IgG     2   1.000000    0.0
6  IgA,IgG     IgG    26   3.750000    0.0
7  IgD,IgG     IgG     8   1.400000    0.0
8      IgG     IgA     2   5.200000    1.0
9      IgG IgD,IgG     1   1.000000    0.0
10     IgG     IgG   136 162.500000    1.0

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.



