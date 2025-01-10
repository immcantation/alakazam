**makeChangeoClone** - *Generate a ChangeoClone object for lineage construction*

Description
--------------------

`makeChangeoClone` takes a data.frame with AIRR or Change-O style columns as input and 
masks gap positions, masks ragged ends, removes duplicates sequences, and merges 
annotations associated with duplicate sequences. It returns a `ChangeoClone` 
object which serves as input for lineage reconstruction. **Note**: To use the 
most recent methods for building, visualizing and analyzing 
trees, use the R package [Dowser](https://dowser.readthedocs.io).


Usage
--------------------
```
makeChangeoClone(
data,
id = "sequence_id",
seq = "sequence_alignment",
germ = "germline_alignment",
v_call = "v_call",
j_call = "j_call",
junc_len = "junction_length",
clone = "clone_id",
mask_char = "N",
locus = "locus",
max_mask = 0,
pad_end = FALSE,
text_fields = NULL,
num_fields = NULL,
seq_fields = NULL,
add_count = TRUE,
verbose = FALSE
)
```

Arguments
-------------------

data
:   data.frame containing the AIRR or Change-O data for a clone. See Details
for the list of required columns and their default values.

id
:   name of the column containing sequence identifiers.

seq
:   name of the column containing observed DNA sequences. All 
sequences in this column must be multiple aligned.

germ
:   name of the column containing germline DNA sequences. All entries 
in this column should be identical for any given clone, and they
must be multiple aligned with the data in the `seq` column.

v_call
:   name of the column containing V-segment allele assignments. All 
entries in this column should be identical to the gene level.

j_call
:   name of the column containing J-segment allele assignments. All 
entries in this column should be identical to the gene level.

junc_len
:   name of the column containing the length of the junction as a 
numeric value. All entries in this column should be identical 
for any given clone.

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

mask_char
:   character to use for masking and padding.

locus
:   name of the column containing locus specification. Must be present
and only contain the value "IGH", representing heavy chains.

max_mask
:   maximum number of characters to mask at the leading and trailing
sequence ends. If `NULL` then the upper masking bound will 
be automatically determined from the maximum number of observed 
leading or trailing Ns amongst all sequences. If set to `0` 
(default) then masking will not be performed.

pad_end
:   if `TRUE` pad the end of each sequence with `mask_char`
to make every sequence the same length.

text_fields
:   text annotation columns to retain and merge during duplicate removal.

num_fields
:   numeric annotation columns to retain and sum during duplicate removal.

seq_fields
:   sequence annotation columns to retain and collapse during duplicate 
removal. Note, this is distinct from the `seq` and `germ` 
arguments, which contain the primary sequence data for the clone
and should not be repeated in this argument.

add_count
:   if `TRUE` add an additional annotation column called 
`collapse_count` during duplicate removal that indicates the 
number of sequences that were collapsed.

verbose
:   passed on to `collapseDuplicates`. If `TRUE`, report the 
numbers of input, discarded and output sequences; otherwise, process
sequences silently.




Value
-------------------

A [ChangeoClone](ChangeoClone-class.md) object containing the modified clone.


Details
-------------------

The input data.frame (`data`) must columns for each of the required column name 
arguments: `id`, `seq`, `germ`, `v_call`, `j_call`, 
`junc_len`, and `clone`.  The default values are as follows:

+ `id       = "sequence_id"`:         unique sequence identifier.
+ `seq      = "sequence_alignment"`:  IMGT-gapped sample sequence.
+ `germ     = "germline_alignment"`:  IMGT-gapped germline sequence.
+ `v_call   = "v_call"`:              V segment allele call.
+ `j_call   = "j_call"`:              J segment allele call.
+ `junc_len = "junction_length"`:     junction sequence length.
+ `clone    = "clone_id"`:            clone identifier.

Additional annotation columns specified in the `text_fields`, `num_fields` 
or `seq_fields` arguments will be retained in the `data` slot of the return 
object, but are not required. If the input data.frame `data` already contains a 
column named `sequence`, which is not used as the `seq` argument, then that 
column will not be retained.

The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
However, all sequences (both observed and germline) must be multiple aligned using
some scheme for both proper duplicate removal and lineage reconstruction. 

The value for the germline sequence, V-segment gene call, J-segment gene call, 
junction length, and clone identifier are determined from the first entry in the 
`germ`, `v_call`, `j_call`, `junc_len` and `clone` columns, 
respectively. For any given clone, each value in these columns should be identical.



Examples
-------------------

```R
# Example data
db <- data.frame(sequence_id=LETTERS[1:4],
sequence_alignment=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
germline_alignment="CCCCAGGG",
v_call="Homsap IGKV1-39*01 F",
j_call="Homsap IGKJ5*01 F",
junction_length=2,
clone_id=1,
locus=rep("IGH", length=4),
c_call=c("IGHM", "IGHG", "IGHG", "IGHA"),
duplicate_count=1:4,
stringsAsFactors=FALSE)


 # Without end masking
 makeChangeoClone(db, text_fields="c_call", num_fields="duplicate_count")

```


```
An object of class "ChangeoClone"
Slot "data":
  sequence_id sequence    c_call duplicate_count collapse_count
1           C NAACTGGN      IGHG               3              1
2           A CCCCTGGG IGHG,IGHM               3              2

Slot "clone":
[1] "1"

Slot "germline":
[1] "CCCCAGGG"

Slot "v_gene":
[1] "IGKV1-39"

Slot "j_gene":
[1] "IGKJ5"

Slot "junc_len":
[1] 2


```


```R

 # With end masking
 makeChangeoClone(db, max_mask=3, text_fields="c_call", num_fields="duplicate_count")

```


```
An object of class "ChangeoClone"
Slot "data":
  sequence_id sequence         c_call duplicate_count collapse_count
1           A NNNCTGNN IGHA,IGHG,IGHM              10              4

Slot "clone":
[1] "1"

Slot "germline":
[1] "CCCCAGGG"

Slot "v_gene":
[1] "IGKV1-39"

Slot "j_gene":
[1] "IGKJ5"

Slot "junc_len":
[1] 2


```



See also
-------------------

Executes in order [maskSeqGaps](maskSeqGaps.md), [maskSeqEnds](maskSeqEnds.md), 
[padSeqEnds](padSeqEnds.md), and [collapseDuplicates](collapseDuplicates.md). 
Returns a [ChangeoClone](ChangeoClone-class.md) object which serves as input to
[buildPhylipLineage](buildPhylipLineage.md).






