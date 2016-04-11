





**DiversityCurve-class** - *S4 class defining diversity curve*

Description
--------------------

`DiversityCurve` defines diversity (<code class = 'eq'>D</code>) scores over multiple diversity 
orders (<code class = 'eq'>Q</code>).

Usage
--------------------

```
"print"(x)
```

Arguments
-------------------

x
:   DiversityCurve object



Slots
-------------------



`data`
data.frame defining the diversity curve with the following columns:

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



`groups`
character vector of groups retained in the diversity calculation.


`n`
numeric vector indication the number of sequences sampled from each group.


`nboot`
number of bootstrap realizations performed.


`ci`
confidence interval defining the upper and lower bounds 
(a value between 0 and 1).






