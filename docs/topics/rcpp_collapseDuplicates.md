





**rcpp_collapseDuplicates** - *Remove duplicate DNA sequences and combine annotations*

Description
--------------------

`collapseDuplicates` identifies duplicate DNA sequences, allowing for ambiguous 
characters, removes the duplicate entries, and combines any associated annotations.


Usage
--------------------
```
rcpp_collapseDuplicates(data, id = "SEQUENCE_ID", seq = "SEQUENCE_IMGT",
text_fields = NULL, num_fields = NULL, seq_fields = NULL,
add_count = FALSE, ignore = c("N", "-", ".", "?"), sep = ",",
verbose = FALSE)
```

Arguments
-------------------

data
:   data.frame containing Change-O columns. The data.frame 
must contain, at a minimum, a unique identifier column 
and a column containg a character vector of DNA sequences.

id
:   name of the column containing sequence identifiers.

seq
:   name of the column containing DNA sequences.

text_fields
:   character vector of textual columns to collapse. The textual 
annotations of duplicate sequences will be merged into a single 
string with each unique value alphabetized and delimited by 
`sep`.

num_fields
:   vector of numeric columns to collapse. The numeric annotations
of duplicate sequences will be summed.

seq_fields
:   vector of nucletoide sequence columns to collapse. The sequence 
with the fewest numer of non-informative characters will be retained. 
Where a non-informative character is one of `c("N", "-", ".", "?")`. 
Note, this is distinct from the `seq` parameter which is used to 
determine duplicates.

add_count
:   if `TRUE` add the column `COLLAPSE_COUNT` that 
indicates the number of sequences that were collapsed to build 
each unique entry.

ignore
:   vector of characters to ignore when testing for equality.

sep
:   character to use for delimiting collapsed annotations in the 
`text_fields` columns. Defines both the input and output 
delimiter.

verbose
:   if `TRUE` report the number input, discarded and output 
sequences; if `FALSE` process sequences silently.



Value
-------------------

A modified `data` data.frame with duplicate sequences removed and 
annotation fields collapsed.

Details
-------------------

`collapseDuplicates` identifies duplicate sequences in the `seq` column by
testing for character identity, with consideration of IUPAC ambiguous nucleotide codes. 
A cluster of sequences are considered duplicates if they are all equivalent, and no 
member of the cluster is equivalent to a sequence in a different cluster. 

Textual annotations, specified by `text_fields`, are collapsed by taking the unique
set of values within in each duplicate cluster and delimiting those values by `sep`.
Numeric annotations, specified by `num_fields`, are collapsed by summing all values 
in the duplicate cluster. Sequence annotations, specified by `seq_fields`, are 
collapsed by retaining the first sequence with the fewest number of N characters.

Columns that are not specified in either `text_fields`, `num_fields`, or 
`seq_fields` will be retained, but the value will be chosen from a random entry 
amongst all sequences in a cluster of duplicates.

An ambiguous sequence is one that can be assigned to two different clusters, wherein
the ambiguous sequence is equivalent to two sequences which are themselves 
non-equivalent. Ambiguous sequences arise due to ambiguous characters at positions that
vary across sequences, and are discarded along with their annotations. Thus, ambiguous
sequences are removed as duplicates of some sequence, but do not create a potential
false-positive annotation merger. Ambiguous sequences are not included in the 
`COLLAPSE_COUNT` annotation that is added when `add_count=TRUE`.



Examples
-------------------

```R
# Example Change-O data.frame
df <- data.frame(SEQUENCE_ID=LETTERS[1:4],
SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
TYPE=c("IgM", "IgG", "IgG", "IgA"),
SAMPLE=c("S1", "S1", "S2", "S2"),
COUNT=1:4,
stringsAsFactors=FALSE)

# Annotations are not parsed if neither text_fields nor num_fields is specified
# The retained sequence annotations will be random
collapseDuplicates(df, verbose=TRUE)

```


```
 FUNCTION> collapseDuplicates
    TOTAL> 4
   UNIQUE> 2
COLLAPSED> 1
DISCARDED> 1


```


```
  SEQUENCE_ID SEQUENCE_IMGT TYPE SAMPLE COUNT
1           C      NAACTGGN  IgG     S2     3
2           A      CCCCTGGG  IgM     S1     1

```


```R

# Unique text_fields annotations are combined into a single string with ","
# num_fields annotations are summed
# Ambiguous duplicates are discarded
collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
verbose=TRUE)

```


```
 FUNCTION> collapseDuplicates
    TOTAL> 4
   UNIQUE> 2
COLLAPSED> 1
DISCARDED> 1


```


```
  SEQUENCE_ID SEQUENCE_IMGT    TYPE SAMPLE COUNT
1           C      NAACTGGN     IgG     S2     3
2           A      CCCCTGGG IgG,IgM     S1     3

```


```R

# Use alternate delimiter for collapsing textual annotations
collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
sep="/", verbose=TRUE)

```


```
 FUNCTION> collapseDuplicates
    TOTAL> 4
   UNIQUE> 2
COLLAPSED> 1
DISCARDED> 1


```


```
  SEQUENCE_ID SEQUENCE_IMGT    TYPE SAMPLE COUNT
1           C      NAACTGGN     IgG     S2     3
2           A      CCCCTGGG IgG/IgM     S1     3

```


```R

# Add count of duplicates
collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
add_count=TRUE, verbose=TRUE)

```


```
 FUNCTION> collapseDuplicates
    TOTAL> 4
   UNIQUE> 2
COLLAPSED> 1
DISCARDED> 1


```


```
  SEQUENCE_ID SEQUENCE_IMGT    TYPE SAMPLE COUNT COLLAPSE_COUNT
1           C      NAACTGGN     IgG     S2     3              1
2           A      CCCCTGGG IgG,IgM     S1     3              2

```


```R

# Masking ragged ends may impact duplicate removal
df$SEQUENCE_IMGT <- maskSeqEnds(df$SEQUENCE_IMGT)
collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
add_count=TRUE, verbose=TRUE)
```


```
 FUNCTION> collapseDuplicates
    TOTAL> 4
   UNIQUE> 1
COLLAPSED> 3
DISCARDED> 0


```


```
  SEQUENCE_ID SEQUENCE_IMGT        TYPE SAMPLE COUNT COLLAPSE_COUNT
1           A      NNNCTGNN IgA,IgG,IgM  S1,S2    10              4

```



See also
-------------------

Equality is tested with [testSeqEqual](testSeqEqual.md).  For IUPAC ambiguous 
character codes see [IUPAC_DNA](IUPAC_CODES.md).


