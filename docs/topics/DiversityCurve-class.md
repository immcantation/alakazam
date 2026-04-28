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

+ `group`:    group label.
+ `q`:        diversity order.
+ `d`:        mean diversity index over all bootstrap 
realizations.
+ `d_sd`:     standard deviation of the diversity index 
over all bootstrap realizations.
+ `d_lower`:  diversity lower confidence interval bound.
+ `d_upper`:  diversity upper confidence interval bound.
+ `e`:        evenness index calculated as `D` 
divided by `D` at `Q=0`.
+ `e_lower`:  evenness lower confidence interval bound.
+ `e_upper`:  evenness upper confidence interval bound.


`tests`
:   data.frame describing the significance test results with columns:

+ `test`:        string listing the two groups tested.
+ `delta_mean`:  mean of the <code class = 'eq'>D</code> bootstrap delta 
distribution for the test.
+ `delta_sd`:    standard deviation of the <code class = 'eq'>D</code> 
bootstrap delta distribution for the test.
+ `pvalue`:      p-value for the test.


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









