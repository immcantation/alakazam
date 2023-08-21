**getMRCA** - *Retrieve the first non-root node of a lineage tree*

Description
--------------------

`getMRCA` returns the set of lineage tree nodes with the minimum weighted or 
unweighted path length from the root (germline) of the lineage tree, allowing for 
exclusion of specific groups of nodes.


Usage
--------------------
```
getMRCA(
graph,
path = c("distance", "steps"),
root = "Germline",
field = NULL,
exclude = NULL
)
```

Arguments
-------------------

graph
:   igraph object containing an annotated lineage tree.

path
:   string defining whether to use unweighted (steps) or weighted (distance) 
measures for determining the founder node set..

root
:   name of the root (germline) node.

field
:   annotation field to use for both unweighted path length exclusion and
consideration as an MRCA node. If `NULL` do not exclude any nodes.

exclude
:   vector of annotation values in `field` to exclude from the potential 
MRCA set. If `NULL` do not exclude any nodes. Has no effect if 
`field=NULL`.




Value
-------------------

A data.frame of the MRCA node(s) containing the columns:

+ `name`:      node name
+ `steps`:     path length as the number of nodes traversed
+ `distance`:  path length as the sum of edge weights

Along with additional columns corresponding to the 
annotations of the input graph.



Examples
-------------------

```R
# Define example graph
graph <- ExampleTrees[[23]]

# Use unweighted path length and do not exclude any nodes
getMRCA(graph, path="steps", root="Germline")

```

*This graph was created by an old(er) igraph version.
  Call upgrade_graph() on it to use with the current igraph version
  For now we convert it on the fly...*
```
                         name
GN5SHBT06HH3QD GN5SHBT06HH3QD
                                                                                                                                                                                                                                                                                                                                                                                                                       sequence
GN5SHBT06HH3QD GAGGTGCAGCTGGTGGTATCTGGGGGANNNGGCTTGGTACAGCCAGGGCGGTCCCTAAGACTCTCCTGTACAGTTTCTGGATTCACCTTTNNNNNNNNNNNNGGTGATTATGCTATGACCTGGATCCGCCAGGCTCCTGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAGCAAAACTTTTGGTGGGACAGCAGATTACGCCGCGTTTGTGAGANNNGGCAGATTCACCATCTCAAGAGATGATTCCAAAAACATCGCCTATCTGCAATTGAACAGCCTGAAAACCGAGGACACAGGCGTCTATTACTGTGGTAGGGATCTCGCCGTAAGTGACACAATAGGTGGTACTAACTGGTTCGACCCCTGGGGCCAGGGGACCCCGGTCACCGTCTCCTCAG
               sample_id c_call duplicate_count          label steps distance
GN5SHBT06HH3QD       +7d   IGHA              10 GN5SHBT06HH3QD     1       20

```


```R

# Exclude nodes without an isotype annotation and use weighted path length
getMRCA(graph, path="distance", root="Germline", field="c_call", exclude=NA)
```


```
                         name
GN5SHBT06HH3QD GN5SHBT06HH3QD
                                                                                                                                                                                                                                                                                                                                                                                                                       sequence
GN5SHBT06HH3QD GAGGTGCAGCTGGTGGTATCTGGGGGANNNGGCTTGGTACAGCCAGGGCGGTCCCTAAGACTCTCCTGTACAGTTTCTGGATTCACCTTTNNNNNNNNNNNNGGTGATTATGCTATGACCTGGATCCGCCAGGCTCCTGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAGCAAAACTTTTGGTGGGACAGCAGATTACGCCGCGTTTGTGAGANNNGGCAGATTCACCATCTCAAGAGATGATTCCAAAAACATCGCCTATCTGCAATTGAACAGCCTGAAAACCGAGGACACAGGCGTCTATTACTGTGGTAGGGATCTCGCCGTAAGTGACACAATAGGTGGTACTAACTGGTTCGACCCCTGGGGCCAGGGGACCCCGGTCACCGTCTCCTCAG
               sample_id c_call duplicate_count          label steps distance
GN5SHBT06HH3QD       +7d   IGHA              10 GN5SHBT06HH3QD     1       20

```



See also
-------------------

Path lengths are determined with [getPathLengths](getPathLengths.md).






