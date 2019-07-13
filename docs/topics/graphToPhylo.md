**graphToPhylo** - *Convert a tree in igraph `graph` format to ape `phylo` format.*

Description
--------------------

`graphToPhylo` a tree in igraph `graph` format to ape `phylo` 
format.


Usage
--------------------
```
graphToPhylo(graph)
```

Arguments
-------------------

graph
:   An igraph `graph` object.




Value
-------------------

A `phylo` object representing the input tree. Tip and internal node names are 
stored in the `tip.label` and `node.label` vectors, respectively.


Details
-------------------

Convert from igraph `graph` object to ape `phylo` object. If `graph` object
was previously rooted with the germline as the direct ancestor, this will re-attach the 
germline as a descendant node with a zero branch length to a new universal common ancestor (UCA) 
node and store the germline node ID in the `germid` attribute and UCA node number in 
the `uca` attribute. Otherwise these attributes will not be specified in the `phylo` object. 
Using `phyloToGraph(phylo, germline=phylo$germid)` creates a `graph` object with the germline 
back as the direct ancestor. Tip and internal node names are 
stored in the `tip.label` and `node.label` vectors, respectively.


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

```

*
Attaching package: ‘igraph’
**The following object is masked from ‘package:testthat’:

    compare
**The following objects are masked from ‘package:dplyr’:

    as_data_frame, groups, union
**The following objects are masked from ‘package:stats’:

    decompose, spectrum
**The following object is masked from ‘package:base’:

    union
*
```R
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




