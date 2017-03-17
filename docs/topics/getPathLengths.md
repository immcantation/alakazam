





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
# Define example graph
graph <- ExampleTrees[[24]]

# Consider all nodes
getPathLengths(graph, root="Germline")

```


```
            NAME STEPS DISTANCE
1      Inferred1     1       20
2 GN5SHBT04CW57C     2       26
3      Inferred2     3       28
4 GN5SHBT08I7RKL     4       29
5 GN5SHBT04CAVIG     5       30
6       Germline     0        0
7 GN5SHBT01D6X0W     2       22
8 GN5SHBT06H7TQD     6       31
9 GN5SHBT05HEG2J     4       33

```


```R

# Exclude nodes without an isotype annotation from step count
getPathLengths(graph, root="Germline", field="ISOTYPE", exclude=NA)
```


```
            NAME STEPS DISTANCE
1      Inferred1     0       20
2 GN5SHBT04CW57C     1       26
3      Inferred2     1       28
4 GN5SHBT08I7RKL     2       29
5 GN5SHBT04CAVIG     3       30
6       Germline     0        0
7 GN5SHBT01D6X0W     1       22
8 GN5SHBT06H7TQD     4       31
9 GN5SHBT05HEG2J     2       33

```



See also
-------------------

See [buildPhylipLineage](buildPhylipLineage.md) for generating input trees.



