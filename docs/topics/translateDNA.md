





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
library(alakazam)
# Load Change-O file
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

translateDNA(df$JUNCTION[1:3])

```


```
[1] "CARDRSTPWRRGIASTTVRTSW" "CARDLLWSVLLTGYYSYGMDAW" "CARDLLWSVLLTGYYSYGMDAW"

```


```R
translateDNA(df$JUNCTION[1:3], trim=TRUE)

```


```
[1] "ARDRSTPWRRGIASTTVRTS" "ARDLLWSVLLTGYYSYGMDA" "ARDLLWSVLLTGYYSYGMDA"

```


```R
translateDNA("ACTGACTCGA")
```


```
[1] "TDS"

```



See also
-------------------

`[translate](http://www.inside-r.org/packages/cran/seqinr/docs/translate)`.



