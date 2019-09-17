**combineIgphyml** - *Combine IgPhyML object parameters into a dataframe*

Description
--------------------

`combineIgphyml` combines IgPhyML object parameters into a data.frame.


Usage
--------------------
```
combineIgphyml(iglist, format = c("wide", "long"))
```

Arguments
-------------------

iglist
:   list of igphyml objects returned by [readIgphyml](readIgphyml.md). 
Each must have an `ID` column in its `param` attribute, 
which can be added automatically using the `id` option of 
`readIgphyml`.

format
:   string specifying whether each column of the resulting data.frame
should represent a parameter (`wide`) or if 
there should only be three columns; i.e. ID, varable, and value
(`long`).




Value
-------------------

A data.frame containing HLP model parameter estimates for all igphyml objects.
Only parameters shared among all objects will be returned.


Details
-------------------

`combineIgphyml` combines repertoire-wide parameter estimates from mutliple igphyml
objects produced by readIgphyml into a dataframe that can be easily used for plotting and 
other hypothesis testing analyses.

All igphyml objects used must have an "ID" column in their `param` attribute, which
can be added automatically from the `id` flag of `readIgphyml`.


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
# Read in and combine two igphyml runs
# s1 <- readIgphyml("IB+7d_lineages_gy.tsv_igphyml_stats_hlp.tab", id="+7d")
# s2 <- readIgphyml("IB+7d_lineages_gy.tsv_igphyml_stats_hlp.tab", id="s2")
# combineIgphyml(list(s1, s2))
```



See also
-------------------

[readIgphyml](readIgphyml.md)






