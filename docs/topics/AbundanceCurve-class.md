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



`abund`
:   data.frame with relative clonal abundance data and confidence intervals,
containing the following columns:

+ `GROUP`:  group identifier.
+ `CLONE`:  clone identifier.
+ `P`:      relative abundance of the clone.
+ `LOWER`:  lower confidence inverval bound.
+ `UPPER`:  upper confidence interval bound.
+ `RANK`:   the rank of the clone abundance.


`bootstrap`
:   data.frame of bootstrapped clonal distributions.

`group`
:   string specifying the name of the group column.

`groups`
:   vector specifying the names of unique groups in group column.

`clone`
:   string specifying the name of the clone column.

`nboot`
:   numeric specifying the number of bootstrap iterations to use.

`uniform`
:   TRUE/FALSE specifying whether or not bootstraps were calculated under rarefaction.

`max_n`
:   numeric specifying the number of species used for rarefaction if rarefaction is used.

`min_n`
:   numeric specifying the minumim tolerated number of species in calculation.

`ci`
:   confidence interval defining the upper and lower bounds 
(a value between 0 and 1).





