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

+ `name`:             node name.
+ `parent`:           name of the parent node.
+ `outdegree`:        number of edges leading from the node.
+ `size`:             total number of nodes within the subtree rooted 
at the node.
+ `depth`:            the depth of the subtree that is rooted at 
the node.
+ `pathlength`:       the maximum pathlength beneath the node.
+ `outdegree_norm`:   `outdegree` normalized by the total 
number of edges.
+ `size_norm`:        `size` normalized by the largest
subtree size (the germline).
+ `depth_norm`:       `depth` normalized by the largest
subtree depth (the germline).
+ `pathlength_norm`:  `pathlength` normalized by the largest
subtree pathlength (the germline).

An additional column corresponding to the value of `field` is added when
specified.



Examples
-------------------

```R
# Summarize a tree
graph <- ExampleTrees[[23]]
summarizeSubtrees(graph, fields="c_call", root="Germline")
```


```
            name    c_call         parent outdegree size depth pathlength outdegree_norm
1 GN5SHBT06HH3QD      IGHA       Germline         1    6     3         10      0.1666667
2 GN5SHBT08F45HV IGHA,IGHG GN5SHBT06HH3QD         4    5     2          7      0.6666667
3       Germline      <NA>           <NA>         1    7     4         30      0.1666667
4 GN5SHBT06IFV0R      IGHG GN5SHBT08F45HV         0    1     1          0      0.0000000
5 GN5SHBT08I3P11      IGHG GN5SHBT08F45HV         0    1     1          0      0.0000000
6 GN5SHBT01BXJY7      IGHG GN5SHBT08F45HV         0    1     1          0      0.0000000
7 GN5SHBT01EGEU6      IGHA GN5SHBT08F45HV         0    1     1          0      0.0000000
  size_norm depth_norm pathlength_norm
1 0.8571429       0.75       0.3333333
2 0.7142857       0.50       0.2333333
3 1.0000000       1.00       1.0000000
4 0.1428571       0.25       0.0000000
5 0.1428571       0.25       0.0000000
6 0.1428571       0.25       0.0000000
7 0.1428571       0.25       0.0000000

```



See also
-------------------

See [buildPhylipLineage](buildPhylipLineage.md) for generating input trees. 
See [getPathLengths](getPathLengths.md) for calculating path length to nodes.






