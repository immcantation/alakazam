





**aminoAcidProperties** - *Calculates amino acid chemical properties for sequence data*

Description
--------------------

`aminoAcidProperties` calculates amino acid sequence physicochemical properties, including
length, hydrophobicity, bulkiness, polarity, aliphatic index, net charge, acidic residue
content, basic residue content, and aromatic residue content.

Usage
--------------------

```
aminoAcidProperties(data, property = c("length", "gravy", "bulk", "aliphatic",
"polarity", "charge", "basic", "acidic", "aromatic"), seq = "JUNCTION",
nt = FALSE, trim = FALSE, label = NULL, ...)
```

Arguments
-------------------

data
:   `data.frame` containing sequence data.

property
:   vector strings specifying the properties to be calculated. Defaults
to calculating all defined properties.

seq
:   `character` name of the column containing input 
sequences.

nt
:   boolean, TRUE if the sequences (or sequence) are DNA and will be translated.

trim
:   if `TRUE` remove the first and last codon/amino acids from each
sequence before calculating properties. If `FALSE` do
not modify input sequences.

label
:   name of sequence region to add as prefix to output column names.

...
:   additional named arguments to pass to the functions 
[gravy](gravy.md), [bulk](bulk.md), [aliphatic](aliphatic.md), [polar](polar.md) or [charge](charge.md).



Value
-------------------

A modified `data` data.frame with the following columns:

+ `*_AA_LENGTH`:     number of amino acids.
+ `*_AA_GRAVY`:      grand average of hydrophobicity (GRAVY) index.
+ `*_AA_BULK`:       average bulkiness of amino acids.
+ `*_AA_ALIPHATIC`:  aliphatic index.
+ `*_AA_POLARITY`:   average polarity of amino acids.
+ `*_AA_CHARGE`:     normalized net charge.
+ `*_AA_BASIC`:      fraction of informative positions that are 
Arg, His or Lys.
+ `*_AA_ACIDIC`:     fraction of informative positions that are 
Asp or Glu.
+ `*_AA_AROMATIC`:   fraction of informative positions that are 
His, Phe, Trp or Tyr.



Where `*` is the value from `label` or the name specified for 
`seq` if `label=NULL`.

Details
-------------------

For all properties except for length, non-informative positions are excluded, 
where non-informative is defined as any character in `c("X", "-", ".", "*")`.

The scores for GRAVY, bulkiness and polarity are calculated as simple averages of the 
scores for each informative positions. The basic, acid and aromatic indices are 
calculated as the fraction of informative positions falling into the given category.

The aliphatic index is calculated using the Ikai, 1980 method.

The net charge is calculated using the method of Moore, 1985, excluding the N-terminus and
C-terminus charges, and normalizing by the number of informative positions.  The default 
pH for the calculation is 7.4.

The following data sources were used for the default property scores:

+ hydropathy:  Kyte & Doolittle, 1982.  
+ bulkiness:   Zimmerman et al, 1968. 
+ polarity:    Grantham, 1974.
+ pK:          EMBOSS.


References
-------------------


1. Zimmerman JM, Eliezer N, Simha R. The characterization of amino acid sequences 
in proteins by statistical methods. J Theor Biol 21, 170-201 (1968).
1. Grantham R. Amino acid difference formula to help explain protein evolution. 
Science 185, 862-864 (1974).
1. Ikai AJ. Thermostability and aliphatic index of globular proteins. 
J Biochem 88, 1895-1898 (1980).
1. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character 
of a protein. J Mol Biol 157, 105-32 (1982).
1. Moore DS. Amino acid and peptide net charges: A simple calculational procedure. 
Biochem Educ 13, 10-11 (1985).
1. Wu YC, et al. High-throughput immunoglobulin repertoire analysis distinguishes 
between human IgM memory and switched memory B-cell populations. 
Blood 116, 1070-8 (2010).
1. Wu YC, et al. The relationship between CD27 negative and positive B cell 
populations in human peripheral blood. 
Front Immunol 2, 1-12 (2011).
1. [http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html](http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html)




Examples
-------------------

```R
# Load example data
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)
df <- df[c(1,10,100), c("SEQUENCE_ID", "JUNCTION")]

# Calculate default amino acid properties from amino acid sequences
# Use a custom output column prefix
df$JUNCTION_TRANS <- translateDNA(df$JUNCTION)
aminoAcidProperties(df, seq="JUNCTION_TRANS", label="JUNCTION")

```


```
       SEQUENCE_ID                                                                                      JUNCTION
1   GN5SHBT08J26Q4                            TGTGCGAGAGATCGGAGCACGCCCTGGCGGCGTGGGATCGCTTCTACCACGGTACGGACGTCGTGG
10  GNDG01208IKHPT TGTGCGAGAGTCCCCCTTTTTGTAGTGGTGGTAGCTGCTACTCCGTTCGGGGGCCGTCGAACCACTACTACTACTACGGTATGGACGTCTGGG
100 GN5SHBT04B1QX2                                  TGTGTTAGAATTGTAGACCCCCACAGTGGCCGGAATGTCCTGCATGCTGCAGACCTTTGG
                     JUNCTION_TRANS JUNCTION_AA_LENGTH JUNCTION_AA_GRAVY JUNCTION_AA_BULK JUNCTION_AA_ALIPHATIC JUNCTION_AA_POLARITY
1            CARDRSTPWRRGIASTTVRTSW                 22        -0.9181818         14.46227              0.400000             8.550000
10  CARVPLFVVVVAATPFGGRRTTTTTTVWTSG                 31         0.5580645         15.59935              0.783871             7.690323
100            CVRIVDPHSGRNVLHAADLW                 20         0.0600000         15.47300              1.120000             8.270000
    JUNCTION_AA_CHARGE JUNCTION_AA_BASIC JUNCTION_AA_ACIDIC JUNCTION_AA_AROMATIC
1          0.178485857        0.22727273         0.04545455           0.09090909
10         0.094399633        0.09677419         0.00000000           0.09677419
100        0.007533018        0.20000000         0.10000000           0.15000000

```


```R
# Calculate default amino acid properties from DNA sequences
aminoAcidProperties(df, seq="JUNCTION", nt=TRUE)

```


```
       SEQUENCE_ID                                                                                      JUNCTION
1   GN5SHBT08J26Q4                            TGTGCGAGAGATCGGAGCACGCCCTGGCGGCGTGGGATCGCTTCTACCACGGTACGGACGTCGTGG
10  GNDG01208IKHPT TGTGCGAGAGTCCCCCTTTTTGTAGTGGTGGTAGCTGCTACTCCGTTCGGGGGCCGTCGAACCACTACTACTACTACGGTATGGACGTCTGGG
100 GN5SHBT04B1QX2                                  TGTGTTAGAATTGTAGACCCCCACAGTGGCCGGAATGTCCTGCATGCTGCAGACCTTTGG
                     JUNCTION_TRANS JUNCTION_AA_LENGTH JUNCTION_AA_GRAVY JUNCTION_AA_BULK JUNCTION_AA_ALIPHATIC JUNCTION_AA_POLARITY
1            CARDRSTPWRRGIASTTVRTSW                 22        -0.9181818         14.46227              0.400000             8.550000
10  CARVPLFVVVVAATPFGGRRTTTTTTVWTSG                 31         0.5580645         15.59935              0.783871             7.690323
100            CVRIVDPHSGRNVLHAADLW                 20         0.0600000         15.47300              1.120000             8.270000
    JUNCTION_AA_CHARGE JUNCTION_AA_BASIC JUNCTION_AA_ACIDIC JUNCTION_AA_AROMATIC
1          0.178485857        0.22727273         0.04545455           0.09090909
10         0.094399633        0.09677419         0.00000000           0.09677419
100        0.007533018        0.20000000         0.10000000           0.15000000

```


```R

# Use the Grantham, 1974 side chain volume scores from the seqinr package
# Set pH=7.0 for the charge calculation
# Calculate only average volume and charge
# Remove the head and tail amino acids from the junction, thus making it the CDR3
library(seqinr)

```

*Loading required package: ade4*
```R
data(aaindex)
x <- aaindex[["GRAR740103"]]$I
# Rename the score vector to use single-letter codes
names(x) <- translateStrings(names(x), ABBREV_AA)
# Calculate properties
aminoAcidProperties(df, property=c("bulk", "charge"), seq="JUNCTION", nt=TRUE, 
trim=TRUE, label="CDR3", bulkiness=x, pH=7.0)
```


```
       SEQUENCE_ID                                                                                      JUNCTION
1   GN5SHBT08J26Q4                            TGTGCGAGAGATCGGAGCACGCCCTGGCGGCGTGGGATCGCTTCTACCACGGTACGGACGTCGTGG
10  GNDG01208IKHPT TGTGCGAGAGTCCCCCTTTTTGTAGTGGTGGTAGCTGCTACTCCGTTCGGGGGCCGTCGAACCACTACTACTACTACGGTATGGACGTCTGGG
100 GN5SHBT04B1QX2                                  TGTGTTAGAATTGTAGACCCCCACAGTGGCCGGAATGTCCTGCATGCTGCAGACCTTTGG
                     JUNCTION_TRANS CDR3_AA_BULK CDR3_AA_CHARGE
1            CARDRSTPWRRGIASTTVRTSW     73.82500     0.20003889
10  CARVPLFVVVVAATPFGGRRTTTTTTVWTSG     72.58621     0.10344795
100            CVRIVDPHSGRNVLHAADLW     73.25000     0.02678262

```



See also
-------------------

See [countPatterns](countPatterns.md) for counting the occurance of specific amino acid subsequences.
See [gravy](gravy.md), [bulk](bulk.md), [aliphatic](aliphatic.md), [polar](polar.md) and [charge](charge.md) for functions 
that calculate the included properties individually.



