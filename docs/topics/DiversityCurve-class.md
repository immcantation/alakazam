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

Arguments
-------------------

x
:   DiversityCurve object

y
:   ignored.

...
:   arguments to pass to [plotDiversityCurve](plotDiversityCurve.md).




Slots
-------------------



`div`
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


`test`
:   list containing information concerning the result of testing diversity:

+ `test`:    data.frame of p-values from testing.
+ `summary`: data.frame of diversity value means and standard deviations for plotting.


`div_group`
:   string specifying the name of the group column in diversity calculation.

`div_groups`
:   vector specifying the names of unique groups in group column in diversity calculation.

`method`
:   string specifying the type of diversity calculated.

`q`
:   vector of diversity hill diversity indices used for computing diversity.

`ci`
:   confidence interval defining the upper and lower bounds 
(a value between 0 and 1).





