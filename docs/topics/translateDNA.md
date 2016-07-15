





**translateDNA** - *Translate nucleotide sequences to amino acids*

Description
--------------------

`translateDNA` translates nucleotide sequences to amino acid sequences.


Usage
--------------------
```
translateDNA(seq, trim = FALSE)
```

Arguments
-------------------

seq
:   vector of strings defining DNA sequence(s) to be converted to translated.

trim
:   boolean flag to remove 3 nts from both ends of seq
(converts IMGT junction to CDR3 region).



Value
-------------------

A vector of translated sequence strings.



Examples
-------------------

```R
# Translate a single sequence
translateDNA("ACTGACTCGA")

```


```
[1] "TDS"

```


```R

# Translate a vector of sequences
translateDNA(ExampleDb$JUNCTION[1:3])

```


```
[1] "CARVKRRGWRRNSLWFGESTPSDAHRWFDPW" "CARMVILGSCYSRGCGTPGPGDGETQYW"    "CARVGIDVVVPAAIPGFDYYYGMDVW"     

```


```R

# Remove the first and last codon from the translation
translateDNA(ExampleDb$JUNCTION[1:3], trim=TRUE)
```


```
[1] "ARVKRRGWRRNSLWFGESTPSDAHRWFDP" "ARMVILGSCYSRGCGTPGPGDGETQY"    "ARVGIDVVVPAAIPGFDYYYGMDV"     

```



See also
-------------------

`[translate](http://www.inside-r.org/packages/cran/seqinr/docs/translate)`.



