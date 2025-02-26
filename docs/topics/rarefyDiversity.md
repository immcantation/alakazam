**rarefyDiversity** - *Generate a clonal diversity index curve*

Description
--------------------

`rarefyDiversity` divides a set of clones by a group annotation,
resamples the sequences from each group, and calculates diversity
scores (<code class = 'eq'>D</code>) over an interval of diversity orders (<code class = 'eq'>q</code>).


Usage
--------------------
```
rarefyDiversity(
data,
group,
clone = "CLONE",
copy = NULL,
min_q = 0,
max_q = 4,
step_q = 0.05,
min_n = 30,
max_n = NULL,
ci = 0.95,
nboot = 2000,
uniform = TRUE,
cell_id = "cell_id",
progress = FALSE
)
```

Arguments
-------------------

data
:   data.frame with Change-O style columns containing clonal assignments.

group
:   name of the `data` column containing group identifiers.

clone
:   name of the `data` column containing clone identifiers.

copy
:   name of the `data` column containing copy numbers for each 
sequence. If `copy=NULL` (the default), then clone abundance
is determined by the number of sequences. If a `copy` column
is specified, then clone abundances is determined by the sum of 
copy numbers within each clonal group.

min_q
:   minimum value of <code class = 'eq'>q</code>.

max_q
:   maximum value of <code class = 'eq'>q</code>.

step_q
:   value by which to increment <code class = 'eq'>q</code>.

min_n
:   minimum number of observations to sample.
A group with less observations than the minimum is excluded.

max_n
:   maximum number of observations to sample. If `NULL` then no 
maximum is set.

ci
:   confidence interval to calculate; the value must be between 0 and 1.

nboot
:   number of bootstrap realizations to generate.

uniform
:   if `TRUE` then uniformly resample each group to the same 
number of observations. If `FALSE` then allow each group to
be resampled to its original size or, if specified, `max_size`.

cell_id
:   name of the `data` column containing cell identifiers.

progress
:   if `TRUE` show a progress bar.




Value
-------------------

A [DiversityCurve](DiversityCurve-class.md) object summarizing the diversity scores.


Details
-------------------

Clonal diversity is calculated using the generalized diversity index (Hill numbers) 
proposed by Hill (Hill, 1973). See [calcDiversity](calcDiversity.md) for further details.

Diversity is calculated on the estimated complete clonal abundance distribution.
This distribution is inferred by using the Chao1 estimator to estimate the number
of seen clones, and applying the relative abundance correction and unseen clone
frequency described in Chao et al, 2015.

To generate a smooth curve, <code class = 'eq'>D</code> is calculated for each value of <code class = 'eq'>q</code> from
`min_q` to `max_q` incremented by `step_q`.  When `uniform=TRUE`
variability in total sequence counts across unique values in the `group` column 
is corrected by repeated resampling from the estimated complete clonal distribution to a 
common number of sequences.

The diversity index (<code class = 'eq'>D</code>) for each group is the mean value of over all resampling 
realizations. Confidence intervals are derived using the standard deviation of the 
resampling realizations, as described in Chao et al, 2015.


References
-------------------


1. Hill M. Diversity and evenness: a unifying notation and its consequences. 
Ecology. 1973 54(2):427-32.
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
### Not run:
# Group by sample identifier
# div <- rarefyDiversity(ExampleDb, "sample_id", step_q=1, max_q=10, nboot=100)
# plotDiversityCurve(div, legend_title="Sample")
# 
# # Grouping by isotype rather than sample identifier
# div <- rarefyDiversity(ExampleDb, "c_call", min_n=40, step_q=1, max_q=10, 
# nboot=100)
# plotDiversityCurve(div, legend_title="Isotype")

```



See also
-------------------

[alphaDiversity](alphaDiversity.md)






