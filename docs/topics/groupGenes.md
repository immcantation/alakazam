**groupGenes** - *Group sequences by gene assignment*

Description
--------------------

`groupGenes` will group rows by shared V and J gene assignments, 
and optionally also by junction lengths.
Both VH:VL paired single-cell BCR-seq and unpaired bulk-seq (heavy chain-only)
are supported.
In the case of ambiguous (multiple) gene assignments, the grouping may
be specified to be a union across all ambiguous V and J gene pairs, 
analagous to single-linkage clustering (i.e., allowing for chaining).


Usage
--------------------
```
groupGenes(data, v_call = "V_CALL", j_call = "J_CALL",
junc_len = NULL, first = FALSE, v_call_light = NULL,
j_call_light = NULL, junc_len_light = NULL,
separator_within_seq = ",", separator_between_seq = ";")
```

Arguments
-------------------

data
:   data.frame containing sequence data.

v_call
:   name of the column containing the heavy chain V-segment 
allele calls.

j_call
:   name of the column containing the heavy chain J-segment 
allele calls.

junc_len
:   name of the column containing the heavy chain junction
length. Optional.

first
:   if `TRUE` only the first call of the gene assignments 
is used. if `FALSE` the union of ambiguous gene 
assignments is used to group all sequences with any 
overlapping gene calls.

v_call_light
:   name of the column containing the light chain V-segment
allele calls. Only applicable for single-cell mode.

j_call_light
:   name of the column containing the light chain J-segment
allele calls. Only applicable for single-cell mode.

junc_len_light
:   name of the column containing the light chain junction
length. Optional and only applicable for single-cell mode.

separator_within_seq
:   a single character specifying the separator
between multiple annotations for a single
sequence. Defaults to `","`.

separator_between_seq
:   a single character specifying the separator
between multiple sequences. Defaults to `";"`.




Value
-------------------

Returns a modified data.frame with disjoint union indices 
in a new `VJ_GROUP` column. 

Note that if `junc_len` (and `junc_len_light` for single-cell) 
is/are supplied, the grouping this `VJ_GROUP` will have been based on 
V, J, and L simultaneously despite the column name being `VJ_GROUP`.


Details
-------------------

All rows containing `NA` values in their `v_call`, `j_call`, and, if
specified, `v_call_light`, `j_call_light`, `junc_len`, and
`junc_len_light` columns will be removed. A warning will be issued when a row 
containing an `NA` is removed.

By default, ambiguous gene assignments are assumed to be separated by commas. This
can be specified otherwise via `separator_within_seq`.

By supplying `junc_len` (and `junc_len_light` for single-cell), the call
amounts to a 1-stage partitioning of the sequences/cells based on V annotation,
J annotation, and junction length simultaneously. Without supplying these columns,
the call amounts to the first stage of a 2-stage partitioning, in which
sequences/cells are partitioned in the first stage based on V annotation and 
J annotation, and then in the second stage further split based on junction length.


Expectation for single-cell input
-------------------



For single-cell BCR data with VH:VL pairing, it is assumed that 

+  every row represents a single cell
+  each cell possibly contains multiple heavy and/or light chains
+  multiple chains of the same type, if any, are separated by 
`separator_between_seq` without any space
+  every chain has its own V(D)J annotation, in which ambiguous V(D)J 
annotations, if any, are separated by `separator_within_seq`


An example:

+  A cell/row has 1 heavy chain (V and J annotations in `v_call` and `j_call`) 
and 2 light chains (V and J annotations in `v_call_light` and `j_call_light`).
+  V annotations for the light chains: `Homsap IGKV1-39*01 F,Homsap IGKV1D-39*01 F;Homsap IGLV2-11*01 F`.
+  J annotations for the light chains: `Homsap IGLJ2*01 F;Homsap IGLJ3*01 F`.


Notice that for both the V and the J annotations, there is a `";"` separating 
the two annotations for the two light chains.

It cannot be the case that there are 2 V annotations but only 1 J annotation for 
the two light chains. Both J annotations must be spelled out for each light chain,
separated by `separator_between_seq`, even if the annotated alleles are the same.

This one-to-one annotation-to-chain correspondence for both V and J is explicitly
checked and an error raised if the requirement is not met.



Examples
-------------------

```R
# Group by genes
db <- groupGenes(data=ExampleDb)
```




