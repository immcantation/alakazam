**stoufferMeta** - *Weighted meta-analysis of p-values via Stouffer's method*

Description
--------------------

`stoufferMeta` combines multiple weighted p-values into a meta-analysis p-value
using Stouffer's Z-score method.


Usage
--------------------
```
stoufferMeta(p, w = NULL)
```

Arguments
-------------------

p
:   numeric vector of p-values.

w
:   numeric vector of weights.




Value
-------------------

A named numeric vector with the combined Z-score and p-value in the form
`c(Z, pvalue)`.



Examples
-------------------

```R
# Define p-value and weight vectors
p <- c(0.1, 0.05, 0.3)
w <- c(5, 10, 1)

# Unweighted
stoufferMeta(p)

```


```
         Z     pvalue 
1.99232360 0.02316778 

```


```R

# Weighted
stoufferMeta(p, w)

```


```
         Z     pvalue 
2.08291783 0.01862936 

```








