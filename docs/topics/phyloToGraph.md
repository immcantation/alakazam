**phyloToGraph** - *Convert a tree in ape `phylo` format to igraph `graph` format.*

Description
--------------------

`phyloToGraph` converts a tree in `phylo` format to and 
`graph` format.


Usage
--------------------
```
phyloToGraph(phylo, germline = "Germline")
```

Arguments
-------------------

phylo
:   An ape `phylo` object.

germline
:   If specified, places specified tip sequence as the direct 
ancestor of the tree




Value
-------------------

A `graph` object representing the input tree.


Details
-------------------

Convert from phylo to graph object. Uses the node.label vector to label internal nodes. Nodes 
may rotate but overall topology will remain constant.


References
-------------------


1. Hoehn KB, Lunter G, Pybus OG - A Phylogenetic Codon Substitution Model for Antibody 
Lineages. Genetics 2017 206(1):417-427
https://doi.org/10.1534/genetics.116.196303 
 1. Hoehn KB, Vander Heiden JA, Zhou JQ, Lunter G, Pybus OG, Kleinstein SHK - 
Repertoire-wide phylogenetic models of B cell molecular evolution reveal 
evolutionary signatures of aging and vaccination. bioRxiv 2019  
https://doi.org/10.1101/558825 




Examples
-------------------

```R
### Not run:
library(igraph)
# library(ape)
# 
# #convert to phylo
# phylo = graphToPhylo(graph)
# 
# #plot tree using ape
# plot(phylo,show.node.label=TRUE)
# 
# #store as newick tree
# write.tree(phylo,file="tree.newick")
# 
# #read in tree from newick file
# phylo_r = read.tree("tree.newick")
# 
# #convert to igraph
# graph_r = phyloToGraph(phylo_r,germline="Germline")
# 
# #plot graph - same as before, possibly rotated
# plot(graph_r,layout=layout_as_tree)
```








