





**testMRCA** - *Tests for MRCA annotation enrichment in lineage trees*

Description
--------------------

`testMRCA` performs a permutation test on a set of lineage trees to determine
the significance of an annotation's association with the MRCA position of the lineage
trees.


Usage
--------------------
```
testMRCA(graphs, field, root = "Germline", exclude = c("Germline", NA),
nperm = 200, progress = FALSE)
```

Arguments
-------------------

graphs
:   list of igraph object containing annotated lineage trees.

field
:   string defining the annotation field to test.

root
:   name of the root (germline) node.

exclude
:   vector of strings defining `field` values to exclude from the
set of potential founder annotations.

nperm
:   number of permutations to perform.

progress
:   if `TRUE` show a progress bar.



Value
-------------------

An [MRCATest](MRCATest-class.md) object containing the test results and permutation
realizations.



Examples
-------------------

```R
# Perform MRCA test on isotypes
x <- testMRCA(ExampleTrees, "ISOTYPE", nperm=100)

```


```


```


```R
print(x)
```


```
  ANNOTATION COUNT EXPECTED PVALUE
1        IgA    16 13.30000 0.0000
2    IgA,IgG     1  1.34375 0.3125
3        IgG    32 34.77000 1.0000

```



See also
-------------------

Uses [getMRCA](getMRCA.md) and [getPathLengths](getPathLengths.md). 
See [plotMRCATest](plotMRCATest.md) for plotting the permutation distributions.



