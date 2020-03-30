**AbundanceCurve-class** - *S4 class defining a clonal abundance curve*

Description
--------------------

`AbundanceCurve` defines clonal abundance values.


Usage
--------------------
```
"print"(x)
```
```
"plot"(x, y, ...)
```

Arguments
-------------------

x
:   AbundanceCurve object

y
:   ignored.

...
:   arguments to pass to [plotDiversityCurve](plotDiversityCurve.md).




Slots
-------------------



`abundance`
:   data.frame with relative clonal abundance data and confidence intervals,
containing the following columns:

+ `group`:  group identifier.
+ `clone_id` or `CLONE`:  clone identifier. 
+ `p`:      relative abundance of the clone.
+ `lower`:  lower confidence inverval bound.
+ `upper`:  upper confidence interval bound.
+ `rank`:   the rank of the clone abundance.


`bootstrap`
:   data.frame of bootstrapped clonal distributions.

`clone_by`
:   string specifying the name of the clone column.

`group_by`
:   string specifying the name of the grouping column.

`groups`
:   vector specifying the names of unique groups in group column.

`n`
:   numeric vector indication the number of sequences sampled in each group.

`nboot`
:   numeric specifying the number of bootstrap iterations to use.

`ci`
:   confidence interval defining the upper and lower bounds 
(a value between 0 and 1).









