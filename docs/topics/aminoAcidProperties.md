**aminoAcidProperties** - *Calculates amino acid chemical properties for sequence data*

Description
--------------------

`aminoAcidProperties` calculates amino acid sequence physicochemical properties, including
length, hydrophobicity, bulkiness, polarity, aliphatic index, net charge, acidic residue
content, basic residue content, and aromatic residue content.


Usage
--------------------
```
aminoAcidProperties(
data,
property = c("length", "gravy", "bulk", "aliphatic", "polarity", "charge", "basic",
"acidic", "aromatic"),
seq = "junction",
nt = FALSE,
trim = FALSE,
label = NULL,
...
)
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

+ `*_aa_length`:     number of amino acids.
+ `*_aa_gravy`:      grand average of hydrophobicity (gravy) index.
+ `*_aa_bulk`:       average bulkiness of amino acids.
+ `*_aa_aliphatic`:  aliphatic index.
+ `*_aa_polarity`:   average polarity of amino acids.
+ `*_aa_charge`:     net charge.
+ `*_aa_basic`:      fraction of informative positions that are 
Arg, His or Lys.
+ `*_aa_acidic`:     fraction of informative positions that are 
Asp or Glu.
+ `*_aa_aromatic`:   fraction of informative positions that are 
His, Phe, Trp or Tyr.



Where `*` is the value from `label` or the name specified for 
`seq` if `label=NULL`.


Details
-------------------

For all properties except for length, non-informative positions are excluded, 
where non-informative is defined as any character in `c("X", "-", ".", "*")`.

The scores for gravy, bulkiness and polarity are calculated as simple averages of the 
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
# Subset example data
db <- ExampleDb[c(1,10,100), c("sequence_id", "junction")]

# Calculate default amino acid properties from amino acid sequences
# Use a custom output column prefix
db$junction_aa <- translateDNA(db$junction)
aminoAcidProperties(db, seq="junction_aa", label="junction")

```


```
     sequence_id                                                                                      junction
1 GN5SHBT02D2WUN TGTGCGAGAGTCAAGCGAAGAGGTTGGCGAAGGAACTCACTATGGTTCGGGGAGTCCACACCTAGCGATGCCCACCGATGGTTCGACCCCTGG
2 GN5SHBT08JP7HP             TGTGCGAGAGATCGGTATTATTGTGGTGGTGACTGCTATTCCCCCCTACCCCAGTACTACTACTACGGTATGGACGTCTGG
3 GN5SHBT05HH5SE                                     TGTGCGAGTGCCTGTAGCAGTGGTGGCTGCTACGAGGAGAACTGGCTCGACCCCTGG
                      junction_aa junction_aa_length junction_aa_gravy junction_aa_bulk junction_aa_aliphatic
1 CARVKRRGWRRNSLWFGESTPSDAHRWFDPW                 31        -1.2612903         14.72194             0.2838710
2     CARDRYYCGGDCYSPLPQYYYYGMDVW                 27        -0.7037037         14.50222             0.2888889
3             CASACSSGGCYEENWLDPW                 19        -0.3684211         13.18053             0.3105263
  junction_aa_polarity junction_aa_charge junction_aa_basic junction_aa_acidic junction_aa_aromatic
1             8.687097           4.038916        0.25806452         0.09677419            0.2258065
2             7.874074          -1.233769        0.07407407         0.11111111            0.2962963
3             8.284211          -3.221436        0.00000000         0.15789474            0.1578947

```


```R
# Calculate default amino acid properties from DNA sequences
aminoAcidProperties(db, seq="junction", nt=TRUE)

```


```
     sequence_id                                                                                      junction
1 GN5SHBT02D2WUN TGTGCGAGAGTCAAGCGAAGAGGTTGGCGAAGGAACTCACTATGGTTCGGGGAGTCCACACCTAGCGATGCCCACCGATGGTTCGACCCCTGG
2 GN5SHBT08JP7HP             TGTGCGAGAGATCGGTATTATTGTGGTGGTGACTGCTATTCCCCCCTACCCCAGTACTACTACTACGGTATGGACGTCTGG
3 GN5SHBT05HH5SE                                     TGTGCGAGTGCCTGTAGCAGTGGTGGCTGCTACGAGGAGAACTGGCTCGACCCCTGG
                      junction_aa junction_aa_length junction_aa_gravy junction_aa_bulk junction_aa_aliphatic
1 CARVKRRGWRRNSLWFGESTPSDAHRWFDPW                 31        -1.2612903         14.72194             0.2838710
2     CARDRYYCGGDCYSPLPQYYYYGMDVW                 27        -0.7037037         14.50222             0.2888889
3             CASACSSGGCYEENWLDPW                 19        -0.3684211         13.18053             0.3105263
  junction_aa_polarity junction_aa_charge junction_aa_basic junction_aa_acidic junction_aa_aromatic
1             8.687097           4.038916        0.25806452         0.09677419            0.2258065
2             7.874074          -1.233769        0.07407407         0.11111111            0.2962963
3             8.284211          -3.221436        0.00000000         0.15789474            0.1578947

```


```R

# Use the Grantham, 1974 side chain volume scores from the seqinr package
# Set pH=7.0 for the charge calculation
# Calculate only average volume and charge
# Remove the head and tail amino acids from the junction, thus making it the CDR3
library(seqinr)
data(aaindex)
x <- aaindex[["GRAR740103"]]$I
# Rename the score vector to use single-letter codes
names(x) <- translateStrings(names(x), ABBREV_AA)
# Calculate properties
aminoAcidProperties(db, property=c("bulk", "charge"), seq="junction", nt=TRUE, 
trim=TRUE, label="cdr3", bulkiness=x, pH=7.0)
```


```
     sequence_id                                                                                      junction
1 GN5SHBT02D2WUN TGTGCGAGAGTCAAGCGAAGAGGTTGGCGAAGGAACTCACTATGGTTCGGGGAGTCCACACCTAGCGATGCCCACCGATGGTTCGACCCCTGG
2 GN5SHBT08JP7HP             TGTGCGAGAGATCGGTATTATTGTGGTGGTGACTGCTATTCCCCCCTACCCCAGTACTACTACTACGGTATGGACGTCTGG
3 GN5SHBT05HH5SE                                     TGTGCGAGTGCCTGTAGCAGTGGTGGCTGCTACGAGGAGAACTGGCTCGACCCCTGG
                      junction_aa cdr3_aa_bulk cdr3_aa_charge
1 CARVKRRGWRRNSLWFGESTPSDAHRWFDPW     85.00000       4.242920
2     CARDRYYCGGDCYSPLPQYYYYGMDVW     79.76000      -1.064488
3             CASACSSGGCYEENWLDPW     58.79412      -3.058792

```



See also
-------------------

See [countPatterns](countPatterns.md) for counting the occurance of specific amino acid subsequences.
See [gravy](gravy.md), [bulk](bulk.md), [aliphatic](aliphatic.md), [polar](polar.md) and [charge](charge.md) for functions 
that calculate the included properties individually.






