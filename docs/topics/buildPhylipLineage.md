





**buildPhylipLineage** - *Infer an Ig lineage using PHYLIP*

Description
--------------------

`buildPhylipLineage` reconstructs an Ig lineage via maximum parsimony using the 
dnapars application of the PHYLIP package.


Usage
--------------------
```
buildPhylipLineage(clone, dnapars_exec, rm_temp = FALSE, verbose = FALSE)
```

Arguments
-------------------

clone
:   [ChangeoClone](ChangeoClone-class.md) object containing clone data.

dnapars_exec
:   path to the PHYLIP dnapars executable.

rm_temp
:   if `TRUE` delete the temporary directory after running dnapars;
if `FALSE` keep the temporary directory.

verbose
:   if `FALSE` suppress the output of dnapars; 
if `TRUE` STDOUT and STDERR of dnapars will be passed to 
the console.



Value
-------------------

An igraph `graph` object defining the Ig lineage tree. Each unique input 
sequence in `clone` is a vertex of the tree, with additional vertices being
either the germline (root) sequences or inferred intermediates. The `graph` 
object has the following attributes.

Vertex attributes:

+ `name`:      value in the `SEQUENCE_ID` column of the `data` 
slot of the input `clone` for observed sequences. 
The germline (root) vertex is assigned the name 
"Germline" and inferred intermediates are assigned
names with the format "Inferred1", "Inferred2", ....
+ `sequence`:  value in the `SEQUENCE` column of the `data` 
slot of the input `clone` for observed sequences.
The germline (root) vertex is assigned the sequence
in the `germline` slot of the input `clone`.
The sequence of inferred intermediates are extracted
from the dnapars output.
+ `label`:     same as the `name` attribute.

Additionally, each other column in the `data` slot of the input 
`clone` is added as a vertex attribute with the attribute name set to 
the source column name. For the germline and inferred intermediate vertices,
these additional vertex attributes are all assigned a value of `NA`.

Edge attributes:

+ `weight`:    Hamming distance between the `sequence` attributes
of the two vertices.
+ `label`:     same as the `weight` attribute.

Graph attributes:

+ `clone`:     clone identifier from the `clone` slot of the
input `ChangeoClone`.
+ `v_gene`:    V-segment gene call from the `v_gene` slot of 
the input `ChangeoClone`.
+ `j_gene`:    J-segment gene call from the `j_gene` slot of 
the input `ChangeoClone`.
+ `junc_len`:  junction length (nucleotide count) from the 
`junc_len` slot of the input `ChangeoClone`.


Details
-------------------

`buildPhylipLineage` builds the lineage tree of a set of unique Ig sequences via
maximum parsimony through an external call to the dnapars application of the PHYLIP
package. dnapars is called with default algorithm options, except for the search option, 
which is set to "Rearrange on one best tree". The germline sequence of the clone is used 
for the outgroup. 

Following tree construction using dnapars, the dnapars output is modified to allow
input sequences to appear as internal nodes of the tree. Intermediate sequences 
inferred by dnapars are replaced by children within the tree having a Hamming distance 
of zero from their parent node. The distance calculation allows IUPAC ambiguous 
character matches, where an ambiguous character has distance zero to any character in 
the set of characters it represents. Distance calculation and movement of child nodes 
up the tree is repeated until all parent-child pairs have a distance greater than zero 
between them. The germline sequence (outgroup) is moved to the root of the tree and
excluded from the node replacement processes, which permits the trunk of the tree to be
the only edge with a distance of zero. Edge weights of the resultant tree are assigned 
as the distance between each sequence.

References
-------------------


1. Felsenstein J. PHYLIP - Phylogeny Inference Package (Version 3.2). 
Cladistics. 1989 5:164-166.
1. Stern JNH, Yaari G, Vander Heiden JA, et al. B cells populating the multiple 
sclerosis brain mature in the draining cervical lymph nodes. 
Sci Transl Med. 2014 6(248):248ra107.




Examples
-------------------

```R
### Not run:
# Load example data
# file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
# df <- readChangeoDb(file)
# 
# # Preprocess clone
# clone <- subset(df, CLONE == 164)
# clone <- makeChangeoClone(clone, text_fields=c("SAMPLE", "ISOTYPE"), num_fields="DUPCOUNT")
# 
# # Run PHYLIP and process output
# dnapars_exec <- "~/apps/phylip-3.69/dnapars"
# graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)
# 
# # Plot graph with a tree layout
# library(igraph)
# ly <- layout_as_tree(graph, root="Germline", circular=F, flip.y=T)
# plot(graph, layout=ly)
```



See also
-------------------

Takes as input a [ChangeoClone](ChangeoClone-class.md). 
Temporary directories are created with [makeTempDir](makeTempDir.md).
Distance is calculated using [getSeqDistance](getSeqDistance.md). 
See [igraph](http://www.inside-r.org/packages/cran/igraph/docs/aaa-igraph-package) and [igraph.plotting](http://www.inside-r.org/packages/cran/igraph/docs/plot.common) for working 
with igraph `graph` objects.


