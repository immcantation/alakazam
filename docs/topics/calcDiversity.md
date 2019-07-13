**calcDiversity** - *Calculate the diversity index*

Description
--------------------

`calcDiversity` calculates the clonal diversity index for a vector of diversity 
orders.


Usage
--------------------
```
calcDiversity(p, q)
```

Arguments
-------------------

p
:   numeric vector of clone (species) counts or proportions.

q
:   numeric vector of diversity orders.




Value
-------------------

A vector of diversity scores <code class = 'eq'>D</code> for each <code class = 'eq'>q</code>.


Details
-------------------

This method, proposed by Hill (Hill, 1973), quantifies diversity as a smooth function 
(<code class = 'eq'>D</code>) of a single parameter <code class = 'eq'>q</code>. Special cases of the generalized diversity 
index correspond to the most popular diversity measures in ecology: species richness 
(<code class = 'eq'>q = 0</code>), the exponential of the Shannon-Weiner index (<code class = 'eq'>q</code> approaches <code class = 'eq'>1</code>), the 
inverse of the Simpson index (<code class = 'eq'>q = 2</code>), and the reciprocal abundance of the largest 
clone (<code class = 'eq'>q</code> approaches <code class = 'eq'>+\infty</code>). At <code class = 'eq'>q = 0</code> different clones weight equally, 
regardless of their size. As the parameter <code class = 'eq'>q</code> increase from <code class = 'eq'>0</code> to <code class = 'eq'>+\infty</code> 
the diversity index (<code class = 'eq'>D</code>) depends less on rare clones and more on common (abundant) 
ones, thus encompassing a range of definitions that can be visualized as a single curve. 

Values of <code class = 'eq'>q < 0</code> are valid, but are generally not meaningful. The value of <code class = 'eq'>D</code> 
at <code class = 'eq'>q=1</code> is estimated by <code class = 'eq'>D</code> at <code class = 'eq'>q=0.9999</code>.


References
-------------------


1. Hill M. Diversity and evenness: a unifying notation and its consequences. 
Ecology. 1973 54(2):427-32.




Examples
-------------------

```R
# May define p as clonal member counts
p <- c(1, 1, 3, 10)
q <- c(0, 1, 2)
calcDiversity(p, q)

```


```
[1] 4.000000 2.594272 2.027027

```


```R

# Or proportional abundance
p <- c(1/15, 1/15, 1/5, 2/3)
calcDiversity(p, q)
```


```
[1] 4.000000 2.594272 2.027027

```



See also
-------------------

Used by [alphaDiversity](alphaDiversity.md) and [betaDiversity](betaDiversity.md).



