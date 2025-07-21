**testMRCA** - *Tests for MRCA annotation enrichment in lineage trees*

Description
--------------------

`testMRCA` performs a permutation test on a set of lineage trees to determine
the significance of an annotation's association with the MRCA position of the lineage
trees.


Usage
--------------------
```
testMRCA(
graphs,
field,
root = "Germline",
exclude = c("Germline", NA),
nperm = 200,
progress = FALSE
)
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
# Define example tree set
graphs <- ExampleTrees[1:10]

# Perform MRCA test on isotypes
x <- testMRCA(graphs, "c_call", nperm=10)

```


```


```


```R
print(x)

```


```
  annotation count expected pvalue
1       IGHA     4      3.4      0
2       IGHG     6      6.0      0

```



See also
-------------------

Uses [getMRCA](getMRCA.md) and [getPathLengths](getPathLengths.md). 
See [plotMRCATest](plotMRCATest.md) for plotting the permutation distributions.






