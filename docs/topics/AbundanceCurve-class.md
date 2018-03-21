





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



`data`
:   data.frame with relative clonal abundance data and confidence intervals,
containing the following columns:

+ `GROUP`:  group identifier.
+ `CLONE`:  clone identifier.
+ `P`:      relative abundance of the clone.
+ `LOWER`:  lower confidence inverval bound.
+ `UPPER`:  upper confidence interval bound.
+ `RANK`:   the rank of the clone abundance.


`groups`
:   character vector of group values.

`n`
:   numeric vector indication the number of sequences in each group.

`nboot`
:   number of bootstrap realizations performed.

`ci`
:   confidence interval defining the upper and lower bounds 
(a value between 0 and 1).





