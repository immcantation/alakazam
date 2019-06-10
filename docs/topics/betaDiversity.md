**betaDiversity** - *Calculates the pairwise beta diversity*

Description
--------------------

`betaDiversity` takes in a data.frame or [AbundanceCurve](AbundanceCurve-class.md) and computes 
SOMETHING.
Statistical testing for significant pairwise differences in diversity 
is also not currently implemented.


Usage
--------------------
```
betaDiversity(data, comparisons, min_q = 0, max_q = 4, step_q = 0.1,
ci = 0.95, ...)
```

Arguments
-------------------

data
:   data.frame with Change-O style columns containing clonal assignments or
a [AbundanceCurve](AbundanceCurve-class.md) generate by [estimateAbundance](estimateAbundance.md) object 
containing a previously calculated bootstrap distributions of clonal abundance.

comparisons
:   list of comparisons between group members for computing beta diversity.

min_q
:   minimum value of <code class = 'eq'>q</code>.

max_q
:   maximum value of <code class = 'eq'>q</code>.

step_q
:   value by which to increment <code class = 'eq'>q</code>.

ci
:   confidence interval to calculate; the value must be between 0 and 1.

...
:   additional arguments to pass to [estimateAbundance](estimateAbundance.md). Additional arguments
are ignored if a [AbundanceCurve](AbundanceCurve-class.md) is provided as input.




Value
-------------------

A [DiversityCurve](DiversityCurve-class.md) object summarizing the diversity scores.


Details
-------------------

TODO



Examples
-------------------

```R
# TODO
```




