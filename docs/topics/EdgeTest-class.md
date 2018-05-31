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

+ `PARENT`:    parent node annotation.
+ `CHILD`:     child node annotation
+ `COUNT`:     count of observed edges with the given 
parent-child annotation set.
+ `EXPECTED`:  mean count of expected edges for the 
given parent-child relationship.
+ `PVALUE`:    one-sided p-value for the hypothesis that 
the observed edge abundance is greater 
than expected.


`permutations`
:   data.frame containing the raw permutation test data with columns:

+ `PARENT`:  parent node annotation.
+ `CHILD`:   child node annotation
+ `COUNT`:   count of edges with the given parent-child 
annotation set.
+ `ITER`:    numerical index define which permutation
realization each observation corresponds 
to.


`nperm`
:   number of permutation realizations.





