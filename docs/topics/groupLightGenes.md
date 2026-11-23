**groupLightGenes** - *Group sequences by light chain gene assignment*

Description
--------------------

`groupLightGenes` groups sequences within existing clones by shared light 
chain V and J gene assignments, and optionally also by junction lengths. This 
function is intended to be used after clonal assignment to further subdivide 
clones using light chain information. IGH:IGK/IGL and TRB:TRA and TRD:TRG paired 
single-cell BCR/TCR sequencing is supported. Heavy chain sequences lacking a 
paired light chain are assigned to a light chain vj_group based on sequence 
similarity to paired heavy chains within the same clone.


Usage
--------------------
```
groupLightGenes(
data,
v_call = "v_call",
j_call = "j_call",
junc_len = NULL,
sequence_alignment = "sequence_alignment",
cell_id = "cell_id",
clone = "clone_id",
locus = "locus",
first = FALSE,
nproc = 1
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
Used to assign unpaired heavy chains to a light chain
group. Heavy chain sequences with a clone id must not
have `NA` or empty values in this column; the
function stops with an error if they do.

cell_id
:   name of the column containing cell identifiers or barcodes.
If specified, grouping will be performed in single-cell mode
with the behavior governed by the `locus` and
`only_heavy` arguments. If set to `NULL` then the
bulk sequencing data is assumed.

clone
:   The clone id of the dataset as defined by heavy chains.

locus
:   name of the column containing locus information.
Only applicable to single-cell data.
Ignored if `cell_id=NULL`.

first
:   if `TRUE` only the first call of the gene assignments
is used. if `FALSE` the union of ambiguous gene
assignments is used to group all sequences with any
overlapping gene calls.

nproc
:   Number of cores to use. Clones are processed in
parallel with `parallel::mclapply`. Must be a
positive integer. Values above 1 are not supported on
Windows and fall back to 1 with a warning. Warnings
raised while processing a clone are reported once
per distinct message, naming the clones affected,
and an error in any clone stops the function.




Value
-------------------

Returns a modified data.frame with light chain subgroup indices in 
a new `vj_group_light` column. These labels are only unique within
a clone. A combination of clone id and vj_group_light provides true 
clonal subgroups. Heavy chain sequences in clones with no light chain 
data  will have `NA` in `vj_group_light`. Unpaired heavy 
chain sequences within clones that do have light chain data will be 
assigned to the nearest light chain subgroup based on Hamming 
distance to paired heavy chains.

If `junc_len` is supplied, grouping will be based on V gene, 
J gene, and junction length simultaneously. The output column name 
will remain `vj_group_light` regardless.

The output `v_call`, `j_call`, `cell_id`, and `locus`
columns will be converted to type `character` if they were of
type `factor` in the input `data`.


Details
-------------------

Values in the `locus` column must be one of `c("IGH", "IGK", "IGL")` 
for BCR or `c("TRA", "TRB", "TRD", "TRG")` for TCR sequences. Otherwise, 
the function returns an error message and stops.

Heavy chain sequences lacking a paired light chain (i.e., no light chain row shares
the same `cell_id`) are assigned to a vj_group_light value based on 
Hamming distance to paired heavy chain sequences within the same clone. In the 
case of ties, the largest vj_group_light is preferred; if vj_group_light sizes 
are equal, the lowest-numbered vj_group_light is used.

If light chain sequences have NA or "" values for v_call, j_call, or junc_len 
they will not go through the grouping process. Rather, their paired heavy chain
will be treated as if it were a unpaired cell. The resulting grouping will be 
applied to the light chain sequence. 

Cells with ambiguous V or J calls or cells with multiple light chains will be 
transitively linked within a clone, analogous to single-linking clustering. If 
first = FALSE then this transitive linking only happens when there are mutliple
light chains to a cell that fall into different groups. 

Specifying `junc_len` will force `groupLightGenes` to perform grouping
based on V gene, J gene, and junction length simultaneously. If `junc_len=NULL`,
grouping is based on V gene and J gene only.

All rows containing `NA` values in any of the `v_call`, `j_call`,
and `junc_len` (if `junc_len != NULL`) columns will be removed by 
the internal `groupGenes` call. A warning will be issued when this occurs.



Examples
-------------------

```R
# After clonal assignment, further subdivide clones by light chain gene usage
db <- groupLightGenes(ExampleDb, clone = "clone_id", cell_id = "cell_id")

```

*Warning*:cell_id not found or not specified. Heavy chain sequences will have NA for vj_group_light. For full paired analysis, please specify cell_id.
```R
head(db$vj_group_light)

```


```
[1] NA NA NA NA NA NA

```








