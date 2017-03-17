





**makeChangeoClone** - *Generate a ChangeoClone object for lineage construction*

Description
--------------------

`makeChangeoClone` takes a data.frame with Change-O style columns as input and 
masks gap positions, masks ragged ends, removes duplicates sequences, and merges 
annotations associated with duplicate sequences. It returns a `ChangeoClone` 
object which serves as input for lineage reconstruction.


Usage
--------------------
```
makeChangeoClone(data, id = "SEQUENCE_ID", seq = "SEQUENCE_IMGT",
germ = "GERMLINE_IMGT_D_MASK", vcall = "V_CALL", jcall = "J_CALL",
junc_len = "JUNCTION_LENGTH", clone = "CLONE", max_mask = 0,
text_fields = NULL, num_fields = NULL, seq_fields = NULL,
add_count = TRUE)
```

Arguments
-------------------

data
:   data.frame containing the Change-O data for a clone. See Details
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

vcall
:   name of the column containing V-segment allele assignments. All 
entries in this column should be identical to the gene level.

jcall
:   name of the column containing J-segment allele assignments. All 
entries in this column should be identical to the gene level.

junc_len
:   name of the column containing the length of the junction as a 
numeric value. All entries in this column should be identical 
for any given clone.

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

max_mask
:   maximum number of characters to mask at the leading and trailing
sequence ends. If `NULL` then the upper masking bound will 
be automatically determined from the maximum number of observed 
leading or trailing Ns amongst all sequences. If set to `0` 
(default) then masking will not be performed.

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
COLLAPSE_COUNT during duplicate removal that indicates the 
number of sequences that were collapsed.




Value
-------------------

A [ChangeoClone](ChangeoClone-class.md) object containing the modified clone.


Details
-------------------

The input data.frame (`data`) must columns for each of the required column name 
arguments: `id`, `seq`, `germ`, `vcall`, `jcall`, 
`junc_len`, and `clone`.  The default values are as follows:

+ `id       = "SEQUENCE_ID"`:           unique sequence identifier.
+ `seq      = "SEQUENCE_IMGT"`:         IMGT-gapped sample sequence.
+ `germ     = "GERMLINE_IMGT_D_MASK"`:  IMGT-gapped germline sequence.
+ `vcall    = "V_CALL"`:                V-segment allele call.
+ `jcall    = "J_CALL"`:                J-segment allele call.
+ `junc_len = "JUNCTION_LENGTH"`:       junction sequence length.
+ `clone    = "CLONE"`:                 clone identifier.

Additional annotation columns specified in the `text_fields`, `num_fields` 
or `seq_fields` arguments will be retained in the `data` slot of the return 
object, but are not required. If the input data.frame `data` already contains a 
column named `SEQUENCE`, which is not used as the `seq` argument, then that 
column will not be retained.

The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
However, all sequences (both observed and germline) must be multiple aligned using
some scheme for both proper duplicate removal and lineage reconstruction. 

The value for the germline sequence, V-segment gene call, J-segment gene call, 
junction length, and clone identifier are determined from the first entry in the 
`germ`, `vcall`, `jcall`, `junc_len` and `clone` columns, 
respectively. For any given clone, each value in these columns should be identical.



Examples
-------------------

```R
# Example Change-O data.frame
db <- data.frame(SEQUENCE_ID=LETTERS[1:4],
SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
V_CALL="Homsap IGKV1-39*01 F",
J_CALL="Homsap IGKJ5*01 F",
JUNCTION_LENGTH=2,
GERMLINE_IMGT_D_MASK="CCCCAGGG",
CLONE=1,
TYPE=c("IgM", "IgG", "IgG", "IgA"),
COUNT=1:4,
stringsAsFactors=FALSE)

# Without end masking
makeChangeoClone(db, text_fields="TYPE", num_fields="COUNT")

```


```
An object of class "ChangeoClone"
Slot "data":
  SEQUENCE_ID SEQUENCE    TYPE COUNT COLLAPSE_COUNT
1           C NAACTGGN     IgG     3              1
2           A CCCCTGGG IgG,IgM     3              2

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
makeChangeoClone(db, max_mask=3, text_fields="TYPE", num_fields="COUNT")
```


```
An object of class "ChangeoClone"
Slot "data":
  SEQUENCE_ID SEQUENCE        TYPE COUNT COLLAPSE_COUNT
1           A NNNCTGNN IgA,IgG,IgM    10              4

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

Executes in order [maskSeqGaps](maskSeqGaps.md), [maskSeqEnds](maskSeqEnds.md)
and [collapseDuplicates](collapseDuplicates.md). 
Returns a [ChangeoClone](ChangeoClone-class.md) object which serves as input to
[buildPhylipLineage](buildPhylipLineage.md).



