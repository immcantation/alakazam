





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

+ `ANNOTATION`:  annotation value.
+ `COUNT`:       observed count of MRCA positions 
with the given annotation.
+ `EXPECTED`:    expected mean count of MRCA occurance
for the annotation.
+ `PVALUE`:      one-sided p-value for the hypothesis that 
the observed annotation abundance is greater 
than expected.


`permutations`
:   data.frame containing the raw permutation test data with columns:

+ `ANNOTATION`:  annotation value.
+ `COUNT`:       count of MRCA positions with the 
given annotation.
+ `ITER`:        numerical index define which 
permutation realization each 
observation corresponds to.


`nperm`
:   number of permutation realizations.





