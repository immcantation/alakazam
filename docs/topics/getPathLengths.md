





**getPathLengths** - *Calculate path lengths from the tree root*

Description
--------------------

`getPathLengths` calculates the unweighted (number of steps) and weighted (distance) 
path lengths from the root of a lineage tree.


Usage
--------------------
```
getPathLengths(graph, root = "Germline", field = NULL, exclude = NULL)
```

Arguments
-------------------

graph
:   igraph object containing an annotated lineage tree.

root
:   name of the root (germline) node.

field
:   annotation field to use for exclusion of nodes from step count.

exclude
:   annotation values specifying which nodes to exclude from step count. 
If `NULL` consider all nodes. This does not affect the weighted
(distance) path length calculation.



Value
-------------------

A data.frame with columns:

+ `NAME`:      node name
+ `STEPS`:     path length as the number of nodes traversed
+ `DISTANCE`:  path length as the sum of edge weights




Examples
-------------------

```R
# Define and plot example graph
library(igraph)
graph <- ExampleTrees[[23]]
plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE)

```

![2](getPathLengths-2.png)

```R

# Consider all nodes
getPathLengths(graph, root="Germline")

```


```
            NAME STEPS DISTANCE
1 GN5SHBT06HH3QD     1       20
2 GN5SHBT08F45HV     2       23
3       Germline     0        0
4 GN5SHBT06IFV0R     3       25
5 GN5SHBT08I3P11     3       30
6 GN5SHBT01BXJY7     3       24
7 GN5SHBT01EGEU6     3       24

```


```R

# Exclude nodes without an isotype annotation from step count
getPathLengths(graph, root="Germline", field="ISOTYPE", exclude=NA)
```


```
            NAME STEPS DISTANCE
1 GN5SHBT06HH3QD     1       20
2 GN5SHBT08F45HV     2       23
3       Germline     0        0
4 GN5SHBT06IFV0R     3       25
5 GN5SHBT08I3P11     3       30
6 GN5SHBT01BXJY7     3       24
7 GN5SHBT01EGEU6     3       24

```



See also
-------------------

See [buildPhylipLineage](buildPhylipLineage.md) for generating input trees.



