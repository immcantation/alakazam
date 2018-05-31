**ABBREV_AA** - *Amino acid abbreviation translations*

Description
--------------------

Mappings of amino acid abbreviations.


Usage
--------------------
```
ABBREV_AA
```



Format
-------------------
Named character vector defining single-letter character codes to 
three-letter abbreviation mappings.


Examples
-------------------

```R
aa <- c("Ala", "Ile", "Trp")
translateStrings(aa, ABBREV_AA)
```


```
[1] "A" "I" "W"

```




