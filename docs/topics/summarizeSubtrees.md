





**summarizeSubtrees** - *Generate subtree summary statistics for a tree*

Description
--------------------

`summarizeSubtrees` calculates summary statistics for each node of a tree. Includes
both node properties and subtree properties.


Usage
--------------------
```
summarizeSubtrees(graph, fields = NULL, root = "Germline")
```

Arguments
-------------------

graph
:   igraph object containing an annotated lineage tree.

fields
:   annotation fields to add to the output.

root
:   name of the root (germline) node.



Value
-------------------

A data.frame with columns: 

+ `NAME`:             node name.
+ `PARENT`:           name of the parent node.
+ `OUTDEGREE`:        number of edges leading from the node.
+ `SIZE`:             total number of nodes within the subtree rooted 
at the node.
+ `DEPTH`:            the depth of the subtree that is rooted at 
the node.
+ `PATHLENGTH`:       the maximum pathlength beneath the node.
+ `OUTDEGREE_NORM`:   `OUTDEGREE` normalized by the total 
number of edges.
+ `SIZE_NORM`:        `SIZE` normalized by the largest
subtree size (the germline).
+ `DEPTH_NORM`:       `DEPTH` normalized by the largest
subtree depth (the germline).
+ `PATHLENGTH_NORM`:  `PATHLEGNTH` normalized by the largest
subtree pathlength (the germline).

An additional column corresponding to the value of `field` is added when
specified.



Examples
-------------------

```R
# Define and plot example graph
library(igraph)
graph <- ExampleTrees[[23]]
plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE)

```

![2](summarizeSubtrees-2.png)

```R

# Summarize tree
summarizeSubtrees(graph, fields="ISOTYPE", root="Germline")
```


```
            NAME ISOTYPE         PARENT OUTDEGREE SIZE DEPTH PATHLENGTH OUTDEGREE_NORM SIZE_NORM DEPTH_NORM PATHLENGTH_NORM
1 GN5SHBT06HH3QD     IgA       Germline         1    6     3         10      0.1666667 0.8571429       0.75       0.3333333
2 GN5SHBT08F45HV IgA,IgG GN5SHBT06HH3QD         4    5     2          7      0.6666667 0.7142857       0.50       0.2333333
3       Germline    <NA>           <NA>         1    7     4         30      0.1666667 1.0000000       1.00       1.0000000
4 GN5SHBT06IFV0R     IgG GN5SHBT08F45HV         0    1     1          0      0.0000000 0.1428571       0.25       0.0000000
5 GN5SHBT08I3P11     IgG GN5SHBT08F45HV         0    1     1          0      0.0000000 0.1428571       0.25       0.0000000
6 GN5SHBT01BXJY7     IgG GN5SHBT08F45HV         0    1     1          0      0.0000000 0.1428571       0.25       0.0000000
7 GN5SHBT01EGEU6     IgA GN5SHBT08F45HV         0    1     1          0      0.0000000 0.1428571       0.25       0.0000000

```



See also
-------------------

See [buildPhylipLineage](buildPhylipLineage.md) for generating input trees. 
See [getPathLengths](getPathLengths.md) for calculating path length to nodes.



