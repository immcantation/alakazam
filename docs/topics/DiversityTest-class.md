





**DiversityTest-class** - *S4 class defining diversity significance*

Description
--------------------

`DiversityTest` defines the signifance of diversity (<code class = 'eq'>D</code>) differences at a 
fixed diversity order (<code class = 'eq'>q</code>).


Usage
--------------------
```
"print"(x)
```

Arguments
-------------------

x
:   DiversityTest object



Slots
-------------------



`tests`
:   data.frame describing the significance test results with columns:

+ `test`:          string listing the two groups tested.
+ `pvalue`:        p-value for the test.
+ `delta_mean`:    mean of the <code class = 'eq'>D</code> bootstrap delta 
distribution for the test.
+ `delta_sd`:      standard deviation of the <code class = 'eq'>D</code> 
bootstrap delta distribution for the test.


`summary`
:   data.frame containing summary statistics for the diversity index 
bootstrap distributions, at the given value of <code class = 'eq'>q</code>, with columns:

+ `group`:   the name of the group.
+ `mean`:    mean of the <code class = 'eq'>D</code> bootstrap distribution.
+ `sd`:      standard deviation of the <code class = 'eq'>D</code> bootstrap 
distribution.


`groups`
:   character vector of groups retained in diversity calculation.

`q`
:   diversity order tested (<code class = 'eq'>q</code>).

`n`
:   numeric vector indication the number of sequences sampled from each group.

`nboot`
:   number of bootstrap realizations.





