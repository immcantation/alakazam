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
groupGenes(
data,
v_call = "v_call",
j_call = "j_call",
junc_len = NULL,
cell_id = NULL,
locus = NULL,
only_heavy = TRUE,
first = FALSE
)
```

Arguments
-------------------

data
:   data.frame containing sequence data.

v_call
:   name of the column containing the heavy chain 
V-segment allele calls.

j_call
:   name of the column containing the heavy chain 
J-segment allele calls.

junc_len
:   name of column containing the junction length. Optional.

cell_id
:   name of the column containing cell IDs. Only 
applicable and required for single-cell mode.

locus
:   name of the column containing locus information. 
Only applicable and required for single-cell mode.

only_heavy
:   use only `IGH` (for BCR data) or `TRB/TRD` (for TCR data) 
sequences for grouping. Only applicable and required for 
single-cell mode. Default is `TRUE`.

first
:   if `TRUE` only the first call of the gene assignments 
is used. if `FALSE` the union of ambiguous gene 
assignments is used to group all sequences with any 
overlapping gene calls.




Value
-------------------

Returns a modified data.frame with disjoint union indices 
in a new `vj_group` column. 

Note that if `junc_len` is supplied, the grouping this `vj_group` 
will have been based on V, J, and L simultaneously despite the column name 
being `vj_group`.

Note that the output `v_call`, `j_call`, `cell_id`, and `locus`
columns will be converted to `character` if they were `factor` in the 
input `data`.


Details
-------------------

To invoke single-cell mode, both `cell_id` and `locus` must be supplied. Otherwise,
the function will run under non-single-cell mode, using all input sequences regardless of the
value in the `locus` column.

Values in the `locus` column must be one of `c("IGH", "IGI", "IGK", "IGL"` for BCR 
or `"TRA", "TRB", "TRD", "TRG")` for TCR sequences. Otherwise, the function returns an 
error message and stops.

Under single-cell mode for VH:VL paired sequences, there is a choice of whether grouping
should be done using `IGH` for BCR or `TRB/TRD` for TCR sequences only, or using 
both `IGH, IGK/IGL` for BCR or `TRB/TRD, TRA/TRG` for TCR sequences. 
This is governed by `only_heavy`.

By supplying `junc_len`, the call amounts to a 1-stage partitioning of the sequences/cells 
based on V annotation, J annotation, and junction length simultaneously. Without supplying this 
columns, the call amounts to the first stage of a 2-stage partitioning, in which sequences/cells 
are partitioned in the first stage based on V annotation and J annotation, and then in the second 
stage further split based on junction length.

It is assumed that ambiguous gene assignments are separated by commas.

In the input `data`, the `v_call`, `j_call`, `cell_id`, and `locus` 
columns, if present, must be `character`, as opposed to `factor`.

All rows containing `NA` values in their any of the `v_call`, `j_call`, and, 
if specified, `junc_len`, columns will be removed. A warning will be issued when a row 
containing an `NA` is removed.


Expectation for single-cell input
-------------------



For single-cell BCR data with VH:VL pairing, it is assumed that 

+  every row represents a sequence (chain)
+  heavy and light chains of the same cell are linked by `cell_id`
+  value in `locus` column indicates whether the chain is heavy or light
+  each cell possibly contains multiple heavy and/or light chains
+  every chain has its own V(D)J annotation, in which ambiguous V(D)J 
annotations, if any, are separated by `,` (comma)


An example:

+  A cell has 1 heavy chain and 2 light chains 
+  There should be 3 rows corresponding to this cell
+  One of the light chain has ambiguous V annotation, which looks like `Homsap IGKV1-39*01 F,Homsap IGKV1D-39*01 F`.




Examples
-------------------

```R
# Group by genes
db <- groupGenes(ExampleDb, v_call="v_call", j_call="j_call")
```








