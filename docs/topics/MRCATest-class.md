**MRCATest-class** - *S4 class defining edge significance*

Description
--------------------

`MRCATest` defines the significance of enrichment for annotations appearing at
the MRCA of the tree.


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
:   MRCATest object.

y
:   ignored.

...
:   arguments to pass to [plotMRCATest](plotMRCATest.md).




Slots
-------------------



`tests`
:   data.frame describing the significance test results with columns:

+ `annotation`:  annotation value.
+ `count`:       observed count of MRCA positions 
with the given annotation.
+ `expected`:    expected mean count of MRCA occurrence
for the annotation.
+ `pvalue`:      one-sided p-value for the hypothesis that 
the observed annotation abundance is greater 
than expected.


`permutations`
:   data.frame containing the raw permutation test data with columns:

+ `annotation`:  annotation value.
+ `count`:       count of MRCA positions with the 
given annotation.
+ `iter`:        numerical index define which 
permutation realization each 
observation corresponds to.


`nperm`
:   number of permutation realizations.









