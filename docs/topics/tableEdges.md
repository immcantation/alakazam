**tableEdges** - *Tabulate the number of edges between annotations within a lineage tree*

Description
--------------------

`tableEdges` creates a table of the total number of connections (edges) for each 
unique pair of annotations within a tree over all nodes.


Usage
--------------------
```
tableEdges(graph, field, indirect = FALSE, exclude = NULL)
```

Arguments
-------------------

graph
:   igraph object containing an annotated lineage tree.

field
:   string defining the annotation field to count.

indirect
:   if `FALSE` count direct connections (edges) only. If 
`TRUE` walk through any nodes with annotations specified in 
the `argument` to count indirect connections. Specifying
`indirect=TRUE` with `exclude=NULL` will have no effect.

exclude
:   vector of strings defining `field` values to exclude from counts.
Edges that either start or end with the specified annotations will not
be counted. If `NULL` count all edges.




Value
-------------------

A data.frame defining total annotation connections in the tree with columns:

+ `PARENT`:  parent annotation
+ `CHILD`:   child annotation
+ `COUNT`:   count of edges for the parent-child relationship




Examples
-------------------

```R
# Define example graph
graph <- ExampleTrees[[23]]

# Count direct edges between isotypes including inferred nodes
tableEdges(graph, "ISOTYPE")

```

**Error in matrix(edge_mat, ncol = 2, dimnames = list(NULL, c("PARENT", **: 'data' must be of a vector type, was 'NULL'
```R

# Count direct edges excluding edges to and from germline and inferred nodes
tableEdges(graph, "ISOTYPE", exclude=c("Germline", NA))

```

**Error in matrix(edge_mat, ncol = 2, dimnames = list(NULL, c("PARENT", **: 'data' must be of a vector type, was 'NULL'
```R

# Count indirect edges walking through germline and inferred nodes
tableEdges(graph, "ISOTYPE", indirect=TRUE, exclude=c("Germline", NA))
```

**Error**: Column `PARENT` is unknown

See also
-------------------

See [testEdges](testEdges.md) for performed a permutation test on edge relationships.






