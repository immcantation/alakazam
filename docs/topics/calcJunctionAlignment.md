**calcJunctionAlignment** - *Alignment properties*

Description
--------------------

Report the number of deleted germline nucleotides in the alignment and the
number of V gene and J gene nucleotides in the CDR3.


Usage
--------------------
```
calcJunctionAlignment(
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
:   Rearrangement database

germline_db
:   Reference germline database, for the V,D and J alleles
in `data`

v_call
:   Name of the column containing V-segment allele assignments.

d_call
:   Name of the column containing D-segment allele assignments.

j_call
:   Name of the column containing J-segment allele assignments.

v_germline_start
:   Name of the column containing the number of the starting
position of the alignment in the V reference germline 
in `germline_db`.

v_germline_end
:   Name of the column containing the number of the ending
position of the alignment in the V reference germline 
in `germline_db`.

d_germline_start
:   Name of the column containing the number of the starting
position of the alignment in the D reference germline 
in `germline_db`.

d_germline_end
:   Name of the column containing the number of the ending
position of the alignment in the D reference germline 
in `germline_db`.

j_germline_start
:   Name of the column containing the number of the starting
position of the alignment in the J reference germline 
in `germline_db`.

j_germline_end
:   Name of the column containing the number of the ending
position of the alignment in the J reference germline 
in `germline_db`.

np1_length
:   Combined length of the N and P regions proximal to the V region.

np2_length
:   Combined length of the N and P regions proximal to the J region.

junction
:   Name of the column containing the junction sequence.

junction_length
:   Name of the column containing the length of the 
junction region in nucleotides.

sequence_alignment
:   Name of the column containing the aligned sequence.




Value
-------------------

Six new columns are added to `data`:

1.  v_germline_deleted_3 Number of 3' V germline nucleotides deleted
1.  d_germline_deleted_5 Number of 5' D germline nucleotides deleted
1.  d_germline_deleted_3 Number of 3' D germline nucleotides deleted
1.  j_germline_deleted_5 Number of 5' J germline nucleotides deleted
1.  v_cdr3_length Number of sequence_alignment V nucleotides in the CDR3
1.  j_cdr3_length Number of sequence_alignment J nucleotides in the CDR3




Examples
-------------------

```R
data(oneseq_db)
germline_db <- list(
"IGHV3-11*05"="CAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTCAAGCCTGGAGGGTCCCTGAGACT
CTCCTGTGCAGCCTCTGGATTCACCTTC............AGTGACTACTACATGAGCTGGATCCGCCAGGCTCCAG
GGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTTACACAAACTACGCAGACTCTGTGAAG
...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGA
CACGGCCGTGTATTACTGTGCGAGAGA",
"IGHD3-10*01"="GTATTACTATGGTTCGGGGAGTTATTATAAC",
"IGHJ5*02"="ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
)
oneseq_db <- calcJunctionAlignment(oneseq_db, germline_db)
```








