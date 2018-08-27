**groupGenes** - *Group sequences by gene assignment*

Description
--------------------

`groupGenes` will group rows by shared V and J gene assignments. 
In the case of ambiguous (multiple) gene assignments, the grouping will
be a union across all ambiguous V and J gene pairs, analagous to 
single-linkage clustering (i.e., allowing for chaining).


Usage
--------------------
```
groupGenes(data, v_call = "V_CALL", j_call = "J_CALL", first = FALSE)
```

Arguments
-------------------

data
:   data.frame containing sequence data.

v_call
:   name of the column containing the V-segment allele calls.

j_call
:   name of the column containing the J-segment allele calls.

first
:   if `TRUE` only the first call of the gene assignments 
is used. if `FALSE` the union of ambiguous gene 
assignments is used to group all sequences with any 
overlapping gene calls.




Value
-------------------

Returns a modified `data` data.frame with union indices 
in the `VJ_GROUP` column.


Details
-------------------

All rows containing `NA` valies in their `v_call` or `j_call` column will be removed. 
A warning will be issued when a row containing an `NA` is removed.

Ambiguous gene assignments are assumed to be separated by commas.



Examples
-------------------

```R
# Group by genes
db <- groupGenes(ExampleDb)
```




