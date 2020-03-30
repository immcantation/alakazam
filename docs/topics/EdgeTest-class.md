**EdgeTest-class** - *S4 class defining edge significance*

Description
--------------------

`EdgeTest` defines the significance of parent-child annotation enrichment.


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
:   EdgeTest object.

y
:   ignored.

...
:   arguments to pass to [plotEdgeTest](plotEdgeTest.md).




Slots
-------------------



`tests`
:   data.frame describing the significance test results with columns:

+ `parent`:    parent node annotation.
+ `child`:     child node annotation
+ `count`:     count of observed edges with the given 
parent-child annotation set.
+ `expected`:  mean count of expected edges for the 
given parent-child relationship.
+ `pvalue`:    one-sided p-value for the hypothesis that 
the observed edge abundance is greater 
than expected.


`permutations`
:   data.frame containing the raw permutation test data with columns:

+ `parent`:  parent node annotation.
+ `child`:   child node annotation
+ `count`:   count of edges with the given parent-child 
annotation set.
+ `iter`:    numerical index define which permutation
realization each observation corresponds 
to.


`nperm`
:   number of permutation realizations.









