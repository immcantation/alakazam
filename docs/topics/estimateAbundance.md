**estimateAbundance** - *Estimates the complete clonal relative abundance distribution*

Description
--------------------

`estimateAbundance` estimates the complete clonal relative abundance distribution 
and confidence intervals on clone sizes using bootstrapping.


Usage
--------------------
```
estimateAbundance(data, group, clone = "CLONE", copy = NULL,
uniform = TRUE, nboot = 200, min_n = 30, max_n = NULL,
ci = 0.95)
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

uniform
:   if `TRUE` then uniformly resample each group to the same 
number of observations. If `FALSE` then allow each group to
be resampled to its original size or, if specified, `max_size`.

nboot
:   number of bootstrap realizations to generate.

min_n
:   minimum number of observations to sample.
A group with less observations than the minimum is excluded.

max_n
:   maximum number of observations to sample. If `NULL` then no 
maximum is set.

ci
:   confidence interval to calculate; the value must be between 0 and 1.




Value
-------------------

A [AbundanceCurve](AbundanceCurve-class.md) object summarizing the abundances.


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



