





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
1      IgA     IgA    39  63.020000 0.9700000
2      IgA IgA,IgG     3   3.494949 0.5454545
3      IgA     IgG     2   5.020202 0.7777778
4  IgA,IgG     IgA    29   5.886792 0.0000000
5  IgA,IgG IgA,IgG     1   2.000000 0.6153846
6  IgA,IgG     IgG    24   3.100000 0.0000000
7  IgD,IgG     IgG     8   1.815789 0.0000000
8      IgG     IgA     1   4.380000 0.9500000
9      IgG IgD,IgG     1   1.000000 0.0000000
10     IgG     IgG   115 138.310000 1.0000000

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.



