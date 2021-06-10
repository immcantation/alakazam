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

+ `parent`:  parent annotation
+ `child`:   child annotation
+ `count`:   count of edges for the parent-child relationship




Examples
-------------------

```R
# Define example graph
graph <- ExampleTrees[[23]]

# Count direct edges between isotypes including inferred nodes
tableEdges(graph, "c_call")

```


```
[38;5;246m# A tibble: 4 x 3[39m
[38;5;246m# Groups:   parent [3][39m
  parent    child     count
  [3m[38;5;246m<chr>[39m[23m     [3m[38;5;246m<chr>[39m[23m     [3m[38;5;246m<int>[39m[23m
[38;5;250m1[39m IGHA      IGHA,IGHG     1
[38;5;250m2[39m IGHA,IGHG IGHA          1
[38;5;250m3[39m IGHA,IGHG IGHG          3
[38;5;250m4[39m [31mNA[39m        IGHA          1

```


```R

# Count direct edges excluding edges to and from germline and inferred nodes
tableEdges(graph, "c_call", exclude=c("Germline", NA))

```


```
[38;5;246m# A tibble: 3 x 3[39m
[38;5;246m# Groups:   parent [2][39m
  parent    child     count
  [3m[38;5;246m<chr>[39m[23m     [3m[38;5;246m<chr>[39m[23m     [3m[38;5;246m<int>[39m[23m
[38;5;250m1[39m IGHA      IGHA,IGHG     1
[38;5;250m2[39m IGHA,IGHG IGHA          1
[38;5;250m3[39m IGHA,IGHG IGHG          3

```


```R

# Count indirect edges walking through germline and inferred nodes
tableEdges(graph, "c_call", indirect=TRUE, exclude=c("Germline", NA))
```


```
[38;5;246m# A tibble: 3 x 3[39m
[38;5;246m# Groups:   parent [2][39m
  parent    child     count
  [3m[38;5;246m<chr>[39m[23m     [3m[38;5;246m<chr>[39m[23m     [3m[38;5;246m<int>[39m[23m
[38;5;250m1[39m IGHA      IGHA,IGHG     1
[38;5;250m2[39m IGHA,IGHG IGHA          1
[38;5;250m3[39m IGHA,IGHG IGHG          3

```



See also
-------------------

See [testEdges](testEdges.md) for performed a permutation test on edge relationships.






