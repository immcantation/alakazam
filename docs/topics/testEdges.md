





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
1      IgA     IgA    39  58.000000    1.0
2      IgA IgA,IgG     3   3.200000    0.5
3      IgA     IgG     2   5.300000    0.9
4  IgA,IgG     IgA    29   9.444444    0.0
5  IgA,IgG IgA,IgG     1   2.666667    1.0
6  IgA,IgG     IgG    24   6.250000    0.0
7  IgD,IgG     IgG     8   1.000000    0.0
8      IgG     IgA     1   4.000000    0.9
9      IgG IgD,IgG     1   1.000000    0.0
10     IgG     IgG   112 134.200000    1.0

```



See also
-------------------

Uses [tableEdges](tableEdges.md) and [permuteLabels](permuteLabels.md). 
See [plotEdgeTest](plotEdgeTest.md) for plotting the permutation distributions.


