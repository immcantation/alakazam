





**tableEdges** - *Tabulate the number of edges between annotations within a lineage tree*

Description
--------------------

`tableEdges` creates a table of the total number of connections (edges) for each 
unique pair of annotations within a tree over all nodes.


Usage
--------------------
```
tableEdges(graph, field, indirect = FALSE, exclude = c("Germline", NA))
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
# Define and plot example graph
library(igraph)
graph <- ExampleTrees[[23]]
plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE, 
vertex.size=50, edge.arrow.mode=0, vertex.color="grey80")

```

![2](tableEdges-2.png)

```R

# Count direct edges between isotypes including inferred nodes
tableEdges(graph, "ISOTYPE", exclude=NULL)

```


```
Source: local data frame [4 x 3]
Groups: PARENT [?]

   PARENT   CHILD COUNT
    <chr>   <chr> <int>
1     IgA IgA,IgG     1
2 IgA,IgG     IgA     1
3 IgA,IgG     IgG     3
4    <NA>     IgA     1

```


```R

# Count direct edges excluding edges to and from inferred nodes
tableEdges(graph, "ISOTYPE", exclude=NA)

```


```
Source: local data frame [3 x 3]
Groups: PARENT [?]

   PARENT   CHILD COUNT
    <chr>   <chr> <int>
1     IgA IgA,IgG     1
2 IgA,IgG     IgA     1
3 IgA,IgG     IgG     3

```


```R

# Count indirect edges walking through inferred nodes
tableEdges(graph, "ISOTYPE", indirect=TRUE, exclude=NA)
```


```
Source: local data frame [3 x 3]
Groups: PARENT [?]

   PARENT   CHILD COUNT
    <chr>   <chr> <int>
1     IgA IgA,IgG     1
2 IgA,IgG     IgA     1
3 IgA,IgG     IgG     3

```



See also
-------------------

See [testEdges](testEdges.md) for performed a permutation test on edge relationships.


