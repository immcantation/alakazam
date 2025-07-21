**groupGenes** - *Group sequences by gene assignment*

Description
--------------------

`groupGenes` will group rows by shared V and J gene assignments, 
and optionally also by junction lengths. IGH:IGK/IGL, TRB:TRA, and TRD:TRG 
paired single-cell BCR/TCR sequencing and unpaired bulk sequencing 
(IGH, TRB, TRD chain only) are supported. In the case of ambiguous (multiple) 
gene assignments, the grouping may be specified to be a union across all 
ambiguous V and J gene pairs, analogous to single-linkage clustering 
(i.e., allowing for chaining).


Usage
--------------------
```
groupGenes(
data,
v_call = "v_call",
j_call = "j_call",
junc_len = NULL,
sequence_alignment = NULL,
cell_id = NULL,
split_light = FALSE,
locus = "locus",
only_heavy = TRUE,
first = FALSE
)
```

Arguments
-------------------

data
:   data.frame containing sequence data.

v_call
:   name of the column containing the heavy/long chain 
V-segment allele calls.

j_call
:   name of the column containing the heavy/long chain 
J-segment allele calls.

junc_len
:   name of column containing the junction length.
If `NULL` then 1-stage partitioning is perform
considering only the V and J genes is performed. 
See Details for further clarification.

sequence_alignment
:   name of the column containing the sequence alignment.

cell_id
:   name of the column containing cell identifiers or barcodes. 
If specified, grouping will be performed in single-cell mode
with the behavior governed by the `locus` and 
`only_heavy` arguments. If set to `NULL` then the 
bulk sequencing data is assumed.

split_light
:   A deprecated parameter. This would split clones by the light chain.
For similar function use dowser::resolveLightChains

locus
:   name of the column containing locus information. 
Only applicable to single-cell data.
Ignored if `cell_id=NULL`.

only_heavy
:   This is deprecated. Only heavy chains will be used in clustering.
Use only the IGH (BCR) or TRB/TRD (TCR) sequences 
for grouping. Only applicable to single-cell data.
Ignored if `cell_id=NULL`.

first
:   if `TRUE` only the first call of the gene assignments 
is used. if `FALSE` the union of ambiguous gene 
assignments is used to group all sequences with any 
overlapping gene calls.




Value
-------------------

Returns a modified data.frame with disjoint union indices 
in a new `vj_group` column. 

If `junc_len` is supplied, the grouping this `vj_group` 
will have been based on V, J, and junction length simultaneously. However, 
the output column name will remain `vj_group`.

The output `v_call`, `j_call`, `cell_id`, and `locus`
columns will be converted to type `character` if they were of type 
`factor` in the input `data`.


Details
-------------------

To invoke single-cell mode the `cell_id` argument must be specified and the `locus` 
column must be correct. Otherwise, `groupGenes` will be run with bulk sequencing assumptions, 
using all input sequences regardless of the values in the `locus` column.

Values in the `locus` column must be one of `c("IGH", "IGI", "IGK", "IGL")` for BCR 
or `c("TRA", "TRB", "TRD", "TRG")` for TCR sequences. Otherwise, the function returns an 
error message and stops.

Under single-cell mode with paired chained sequences, there is a choice of whether 
grouping should be done by (a) using IGH (BCR) or TRB/TRD (TCR) sequences only or
(b) using IGH plus IGK/IGL (BCR) or TRB/TRD plus TRA/TRG (TCR). 
This is governed by the `only_heavy` argument.

Specifying `junc_len` will force `groupGenes` to perform a 1-stage partitioning of the 
sequences/cells based on V gene, J gene, and junction length simultaneously. 
If `junc_len=NULL` (no column specified), then `groupGenes` performs only the first 
stage of a 2-stage partitioning in which sequences/cells are partitioned in the first stage 
based on V gene and J gene, and then in the second stage further splits the groups based on 
junction length (the second stage must be performed independently, as this only returns the
first stage results).

In the input `data`, the `v_call`, `j_call`, `cell_id`, and `locus` 
columns, if present, must be of type `character` (as opposed to `factor`). 

It is assumed that ambiguous gene assignments are separated by commas.

All rows containing `NA` values in any of the `v_call`, `j_call`, and `junc_len` 
(if `junc_len != NULL`) columns will be removed. A warning will be issued when a row 
containing an `NA` is removed.


Expectations for single-cell data
-------------------



Single-cell paired chain data assumptions:

+  every row represents a sequence (chain).
+  heavy/long and light/short chains of the same cell are linked by `cell_id`.
+  the value in `locus` column indicates whether the chain is the heavy/long or light/short chain.
+  each cell possibly contains multiple heavy/long and/or light/short chains.
+  every chain has its own V(D)J annotation, in which ambiguous V(D)J 
annotations, if any, are separated by a comma.


Single-cell example:

+  A cell has 1 heavy chain and 2 light chains.
+  There should be 3 rows corresponding to this cell.
+  One of the light chains may have an ambiguous V annotation which looks like `"Homsap IGKV1-39*01 F,Homsap IGKV1D-39*01 F"`.




Examples
-------------------

```R
# Group by genes
db <- groupGenes(ExampleDb)
head(db$vj_group)

```


```
[1] "G28" "G50" "G36" "G36" "G86" "G61"

```








