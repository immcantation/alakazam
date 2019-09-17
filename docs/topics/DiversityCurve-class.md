**DiversityCurve-class** - *S4 class defining a diversity curve*

Description
--------------------

`DiversityCurve` defines diversity (<code class = 'eq'>D</code>) scores over multiple diversity 
orders (<code class = 'eq'>Q</code>).


Usage
--------------------
```
"print"(x)
```
```
"plot"(x, y, ...)
```
```
"plot"(x, y, ...)
```

Arguments
-------------------

x
:   DiversityCurve object

y
:   diversity order to plot (q).

...
:   arguments to pass to [plotDiversityCurve](plotDiversityCurve.md) or [plotDiversityTest](plotDiversityTest.md).




Slots
-------------------



`diversity`
:   data.frame defining the diversity curve with the following columns:

+ `GROUP`:    group label.
+ `Q`:        diversity order.
+ `D`:        mean diversity index over all bootstrap 
realizations.
+ `D_SD`:     standard deviation of the diversity index 
over all bootstrap realizations.
+ `D_LOWER`:  diversity lower confidence inverval bound.
+ `D_UPPER`:  diversity upper confidence interval bound.
+ `E`:        evenness index calculated as `D` 
divided by `D` at `Q=0`.
+ `E_LOWER`:  evenness lower confidence inverval bound.
+ `E_UPPER`:  eveness upper confidence interval bound.


`tests`
:   data.frame describing the significance test results with columns:

+ `TEST`:        string listing the two groups tested.
+ `DELTA_MEAN`:  mean of the <code class = 'eq'>D</code> bootstrap delta 
distribution for the test.
+ `DELTA_SD`:    standard deviation of the <code class = 'eq'>D</code> 
bootstrap delta distribution for the test.
+ `PVALUE`:      p-value for the test.


`group_by`
:   string specifying the name of the grouping column in diversity calculation.

`groups`
:   vector specifying the names of unique groups in group column in diversity calculation.

`method`
:   string specifying the type of diversity calculated.

`q`
:   vector of diversity hill diversity indices used for computing diversity.

`n`
:   numeric vector indication the number of sequences sampled in each group.

`ci`
:   confidence interval defining the upper and lower bounds 
(a value between 0 and 1).









