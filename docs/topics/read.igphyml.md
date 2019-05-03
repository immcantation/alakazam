**read.igphyml** - *Read in output from IgPhyML*

Description
--------------------

`read.igphyml` reads output from the IgPhyML phylogenetics inference package for 
B cell repertoires


Usage
--------------------
```
read.igphyml(file, id = NULL, igraph = TRUE, collapse = TRUE)
```

Arguments
-------------------

file
:   IgPhyML output file (.tab)

id
:   ID to assign to output object

igraph
:   if `TRUE` return trees as igraph `graph` objects. Otherwise,
return trees as ape `phylo` objects.

collapse
:   if `TRUE` transform branch lengths to units of substitutions, 
rather than substitutions per site, and collapse internal nodes
separated by branches < 0.1 substitutions




Value
-------------------

A list containing IgPhyML model parameters and estimated lineage trees. 

Object attributes:

+ `param`:     Dataframe of parameter estimates for each clonal 
lineage. Columns include: `CLONE`, which is the 
clone id; `NSEQ`, the total number of sequences in 
the lineage; `NSITE`, the number of codon sites;
`TREE_LENGTH`, the sum of all branch lengths in 
the estimated lineage tree; and `LHOOD`, the log 
likelihood of the clone's sequences given the tree and
parameters. Subsequent columns are parameter estimates 
from IgPhyML, which will depend on the model used. 
Parameter columns ending with `_MLE` are maximum 
likelihood estimates; those ending with `_LCI` are 
the lower 95
with `_UCI` are the upper 95
estimate. The first line of `param` is for clone 
`REPERTOIRE`, 
which is a summary of all lineages within the repertoire.
For this row, `NSEQ` is the total number of sequences, 
`NSITE` is the average number of sites, and
`TREE_LENGTH` is the mean tree length. For most 
applications, parameter values will be the same for all 
lineages within the repertoire, so access them simply by:
`<object>$parm$OMEGA_CDR_MLE[1]` to, for instance,
get the estimate of dN/dS at the repertoire level.
+ `trees`:     List of tree objects estimated by IgPhyML. If 
`igraph=TRUE` these are igraph `graph` objects. 
If `igraph=TRUE`, these are ape `phylo` objects.
+ `command`:   Command used to run IgPhyML.
+ `id`:        Optional name given to object.



Details
-------------------

`read.igphyml` reads output from the IgPhyML repertoire phylogenetics inference package. 
The resulting object is divded between parameter estimates (usually under the HLP19 model),
which provide information about mutation and selection pressure operating on the sequences.

Trees returned from this function are either igraph objects or phylo objects, and each may be 
visualized accordingly. Futher, branch lengths in tree may represent either the expected number of
substitutions per site (codon, if estimated under HLP or GY94 models), or the total number of 
expected substitutions per site. If the latter, internal nodes - but not tips - separated by branch
lengths less than 0.1 are collapsed to simplify viewing.


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
s1 = read.igphyml("IB+7d_lineages_gy.tsv_igphyml_stats_hlp.tab",id="+7d")

```

*Warning*:cannot open file 'IB+7d_lineages_gy.tsv_igphyml_stats_hlp.tab': No such file or directory**Error in file(file, "rt")**: cannot open the connection
```R
#  print(s1$param$OMEGA_CDR_MLE[1])
#  plot(s1$trees[[1]],layout=layout_as_tree,edge.label=E(s1$trees[[1]])$weight)
```




