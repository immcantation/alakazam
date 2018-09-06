**estimateAbundance** - *Estimates the complete clonal relative abundance distribution*

Description
--------------------

`estimateAbundance` estimates the complete clonal relative abundance distribution 
and confidence intervals on clone sizes using bootstrapping.


Usage
--------------------
```
estimateAbundance(data, group = NULL, clone = "CLONE", copy = NULL,
ci = 0.95, nboot = 2000, progress = FALSE)
```

Arguments
-------------------

data
:   data.frame with Change-O style columns containing clonal assignments.

group
:   name of the `data` column containing group identifiers. 
If `NULL` then no grouping is performed and the `GROUP` 
column of the output will contain the value `NA` for each row.

clone
:   name of the `data` column containing clone identifiers.

copy
:   name of the `data` column containing copy numbers for each 
sequence. If `copy=NULL` (the default), then clone abundance
is determined by the number of sequences. If a `copy` column
is specified, then clone abundances is determined by the sum of 
copy numbers within each clonal group.

ci
:   confidence interval to calculate; the value must be between 0 and 1.

nboot
:   number of bootstrap realizations to generate.

progress
:   if `TRUE` show a progress bar.




Value
-------------------

A data.frame with relative clonal abundance data and confidence intervals,
containing the following columns:

+ `GROUP`:  group identifier. Will be codeNA if `group=NULL`.
+ `CLONE`:  clone identifier.
+ `P`:      relative abundance of the clone.
+ `LOWER`:  lower confidence inverval bound.
+ `UPPER`:  upper confidence interval bound.
+ `RANK`:   the rank of the clone abundance.



Details
-------------------

The complete clonal abundance distribution determined inferred by using the Chao1 
estimator to estimate the number of seen clones, and then applying the relative abundance 
correction and unseen clone frequencies described in Chao et al, 2015.

Confidence intervals are derived using the standard deviation of the resampling 
realizations, as described in Chao et al, 2015.


References
-------------------


1. Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
Scand J Stat. 1984 11, 265270.
1. Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
A framework for sampling and estimation in species diversity studies. 
Ecol Monogr. 2014 84:45-67.
1. Chao A, et al. Unveiling the species-rank abundance distribution by 
generalizing the Good-Turing sample coverage theory. 
Ecology. 2015 96, 11891201.




Examples
-------------------

```R
abund <- estimateAbundance(ExampleDb, "SAMPLE", nboot=100)
```



See also
-------------------

See [plotAbundanceCurve](plotAbundanceCurve.md) for plotting of the abundance distribution.
See [rarefyDiversity](rarefyDiversity.md) for a similar application to clonal diversity.



