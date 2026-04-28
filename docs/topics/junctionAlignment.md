**junctionAlignment** - *Calculate junction region alignment properties*

Description
--------------------

`junctionAlignment` determines the number of deleted germline nucleotides in the 
junction region and the number of V gene and J gene nucleotides in the CDR3.


Usage
--------------------
```
junctionAlignment(
data,
germline_db,
v_call = "v_call",
d_call = "d_call",
j_call = "j_call",
v_germline_start = "v_germline_start",
v_germline_end = "v_germline_end",
d_germline_start = "d_germline_start",
d_germline_end = "d_germline_end",
j_germline_start = "j_germline_start",
j_germline_end = "j_germline_end",
np1_length = "np1_length",
np2_length = "np2_length",
junction = "junction",
junction_length = "junction_length",
sequence_alignment = "sequence_alignment"
)
```

Arguments
-------------------

data
:   `data.frame` containing sequence data.

germline_db
:   reference germline database for the V, D and J genes.
in `data`

v_call
:   V gene assignment column.

d_call
:   D gene assignment column.

j_call
:   J gene assignment column.

v_germline_start
:   column containing the start position of the alignment 
in the V reference germline.

v_germline_end
:   column containing the end position of the alignment in the 
V reference germline.

d_germline_start
:   column containing the start position of the alignment 
in the D reference germline.

d_germline_end
:   column containing the start position of the alignment 
in the D reference germline.

j_germline_start
:   column containing the start position of the alignment 
in the J reference germline.

j_germline_end
:   column containing the start position of the alignment 
in the J reference germline.

np1_length
:   combined length of the N and P regions between the 
V and D regions (heavy chain) or V and J regions (light chain).

np2_length
:   combined length of the N and P regions between the 
D and J regions (heavy chain).

junction
:   column containing the junction sequence.

junction_length
:   column containing the length of the junction region in nucleotides.

sequence_alignment
:   column containing the aligned sequence.




Value
-------------------

A modified input `data.frame` with the following additional columns storing 
junction alignment information:

1. `e3v_length`:     number of 3' V germline nucleotides deleted.
1. `e5d_length`:     number of 5' D germline nucleotides deleted.
1. `e3d_length`:     number of 3' D germline nucleotides deleted.
1. `e5j_length`:     number of 5' J germline nucleotides deleted.
1. `v_cdr3_length`:  number of sequence_alignment V nucleotides in the CDR3.
1. `j_cdr3_length`:  number of sequence_alignment J nucleotides in the CDR3.




Examples
-------------------

```R
germline_db <- list(
"IGHV3-11*05"="CAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTCAAGCCTGGAGGGTCCCTGAGACT
CTCCTGTGCAGCCTCTGGATTCACCTTC............AGTGACTACTACATGAGCTGGATCCGCCAGGCTCCAG
GGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTTACACAAACTACGCAGACTCTGTGAAG
...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGA
CACGGCCGTGTATTACTGTGCGAGAGA",
"IGHD3-10*01"="GTATTACTATGGTTCGGGGAGTTATTATAAC",
"IGHJ5*02"="ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
)

db <- junctionAlignment(SingleDb, germline_db)

```








