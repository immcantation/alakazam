





**testEdges** - *Tests for parent-child annotation enchrichment in lineage trees*

Description
--------------------

`testEdges` performs a permutation test on a set of lineage trees to determine
the significance of an annotation's association with parent-child relationships.


Usage
--------------------
```
testEdges(graphs, field, exclude = c("Germline", NA), nperm = 200,
progress = FALSE)
```

Arguments
-------------------

graphs
:   list of igraph objects with vertex annotations.

field
:   string defining the annotation field to permute.

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
# Perform edge test on isotypes
x <- testEdges(ExampleTrees, "ISOTYPE", nperm=100)

```


```


```


```R
print(x)
```


```
    PARENT   CHILD COUNT   EXPECTED    PVALUE
1      IgA     IgA    39  62.150000 0.9700000
2      IgA IgA,IgG     3   3.353535 0.4848485
3      IgA     IgG     2   4.581633 0.7653061
4  IgA,IgG     IgA    29   6.288136 0.0000000
5  IgA,IgG IgA,IgG     1   1.937500 0.6250000
6  IgA,IgG     IgG    24   3.739130 0.0000000
7  IgD,IgG     IgG     8   1.609756 0.0000000
8      IgG     IgA     1   4.525253 0.9292929
9      IgG IgD,IgG     1   1.000000 0.0000000
10     IgG     IgG   115 138.430000 1.0000000

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.



