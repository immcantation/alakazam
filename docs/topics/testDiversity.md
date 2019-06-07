**testDiversity** - *Pairwise test of the diversity index*

Description
--------------------

`testDiversity` performs pairwise significance tests of the diversity index 
(<code class = 'eq'>D</code>) at a given diversity order (<code class = 'eq'>q</code>) for a set of annotation groups via
rarefaction and bootstrapping.


Usage
--------------------
```
testDiversity(data, q, group, clone = "CLONE", copy = NULL,
min_n = 30, max_n = NULL, nboot = 2000, progress = FALSE,
ci = 0.95)
```

Arguments
-------------------

data
:   data.frame with Change-O style columns containing clonal assignments.

q
:   diversity order to test.

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

min_n
:   minimum number of observations to sample.
A group with less observations than the minimum is excluded.

max_n
:   maximum number of observations to sample. If `NULL` the maximum
if automatically determined from the size of the largest group.

nboot
:   number of bootstrap realizations to perform.

progress
:   if `TRUE` show a progress bar.

ci
:   confidence interval to calculate; the value must be between 0 and 1.




Value
-------------------

A [DiversityCurve](DiversityCurve-class.md) object containing slot test with p-values and summary 
		 statistics.


Details
-------------------

Clonal diversity is calculated using the generalized diversity index proposed by 
Hill (Hill, 1973). See [calcInferredDiversity](calcInferredDiversity.md) for further details.

Diversity is calculated on the estimated complete clonal abundance distribution.
This distribution is inferred by using the Chao1 estimator to estimate the number
of seen clones, and applying the relative abundance correction and unseen clone
frequency described in Chao et al, 2014.

Variability in total sequence counts across unique values in the `group` column is 
corrected by repeated resampling from the estimated complete clonal distribution to 
a common number of sequences. The diversity index estimate (<code class = 'eq'>D</code>) for each group is 
the mean value of over all bootstrap realizations. 

Significance of the difference in diversity index (<code class = 'eq'>D</code>) between groups is tested by 
constructing a bootstrap delta distribution for each pair of unique values in the 
`group` column. The bootstrap delta distribution is built by subtracting the diversity 
index <code class = 'eq'>Da</code> in <code class = 'eq'>group-a</code> from the corresponding value <code class = 'eq'>Db</code> in <code class = 'eq'>group-b</code>, 
for all bootstrap realizations, yeilding a distribution of `nboot` total deltas; where 
<code class = 'eq'>group-a</code> is the group with the greater mean <code class = 'eq'>D</code>. The p-value for hypothesis 
<code class = 'eq'>Da  !=  Db</code> is the value of <code class = 'eq'>P(0)</code> from the empirical cumulative distribution 
function of the bootstrap delta distribution, multiplied by 2 for the two-tailed correction.


Note
-------------------

This method may inflate statistical significance when clone sizes are uniformly small,
such as when most clones sizes are 1, sample size is small, and `max_n` is near
the total count of the smallest data group. Use caution when interpreting the results 
in such cases. We are currently investigating this potential problem.


References
-------------------


1. Hill M. Diversity and evenness: a unifying notation and its consequences. 
Ecology. 1973 54(2):427-32.
1. Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
Scand J Stat. 1984 11, 265270.
1. Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
peripheral blood IgE repertoires in patients with allergic rhinitis. 
J Allergy Clin Immunol. 2014 134(3):604-12.
1. Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
A framework for sampling and estimation in species diversity studies. 
Ecol Monogr. 2014 84:45-67.
1. Chao A, et al. Unveiling the species-rank abundance distribution by 
generalizing the Good-Turing sample coverage theory. 
Ecology. 2015 96, 11891201.




Examples
-------------------

```R
# Groups under the size threshold are excluded and a warning message is issued.
testDiversity(ExampleDb, "SAMPLE", q=0, min_n=30, nboot=100)
```


```
An object of class "DiversityCurve"
Slot "div":
# A tibble: 4 x 9
# Groups:   SAMPLE [2]
  SAMPLE     Q D_ERROR      D D_LOWER D_UPPER     E E_LOWER E_UPPER
  <chr>  <dbl>   <dbl>  <dbl>   <dbl>   <dbl> <dbl>   <dbl>   <dbl>
1 -1h       -1 19058.  39263.  20205.  58321.  49.2  25.3     73.1 
2 -1h        0   279.    798.    518.   1077.   1     0.650    1.35
3 +7d       -1  7726.  10453.   2727.  18179.  54.3  14.2     94.5 
4 +7d        0    67.2   192.    125.    260.   1     0.651    1.35

Slot "test":
$tests
# A tibble: 2 x 5
  test_name  Q     DELTA_MEAN DELTA_SD PVALUE
  <chr>      <chr>      <dbl>    <dbl>  <dbl>
1 -1h != +7d -1        28810.   11142.      0
2 -1h != +7d 0           605.     152.      0

$summary
# A tibble: 4 x 4
# Groups:   SAMPLE [2]
  SAMPLE     Q     SD   MEAN
  <chr>  <dbl>  <dbl>  <dbl>
1 -1h       -1 9723.  39263.
2 -1h        0  143.    798.
3 +7d       -1 3942.  10453.
4 +7d        0   34.3   192.


Slot "method":
[1] "alpha"

Slot "div_group":
[1] "SAMPLE"

Slot "div_groups":
[1] "-1h" "+7d"

Slot "q":
[1] -1  0

Slot "ci":
[1] 0.95


```



See also
-------------------

See [calcInferredDiversity](calcInferredDiversity.md) for the basic calculation.
See [rarefyDiversity](rarefyDiversity.md) for curve generation.
See [ecdf](http://www.rdocumentation.org/packages/stats/topics/ecdf) for computation of the empirical cumulative 
distribution function.



