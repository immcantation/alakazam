**getSegment** - *Get Ig segment allele, gene and family names*

Description
--------------------

`getSegment` performs generic matching of delimited segment calls with a custom 
regular expression. [getAllele](getSegment.md), [getGene](getSegment.md) and [getFamily](getSegment.md) extract 
the allele, gene and family names, respectively, from a character vector of 
immunoglobulin (Ig) or TCR segment allele calls in IMGT format.


Usage
--------------------
```
getSegment(segment_call, segment_regex, first = TRUE, collapse = TRUE,
strip_d = TRUE, omit_nl = FALSE, sep = ",")
```
```
getAllele(segment_call, first = TRUE, collapse = TRUE,
strip_d = TRUE, omit_nl = FALSE, sep = ",")
```
```
getGene(segment_call, first = TRUE, collapse = TRUE, strip_d = TRUE,
omit_nl = FALSE, sep = ",", unique = FALSE)
```
```
getFamily(segment_call, first = TRUE, collapse = TRUE,
strip_d = TRUE, omit_nl = FALSE, sep = ",")
```

Arguments
-------------------

segment_call
:   character vector containing segment calls delimited by commas.

segment_regex
:   string defining the segment match regular expression.

first
:   if `TRUE` return only the first call in 
`segment_call`; if `FALSE` return all calls 
delimited by commas.

collapse
:   if `TRUE` check for duplicates and return only unique 
segment assignments; if `FALSE` return all assignments 
(faster). Has no effect if `first=TRUE`.

strip_d
:   if `TRUE` remove the "D" from the end of gene annotations 
(denoting a duplicate gene in the locus); 
if `FALSE` do not alter gene names.

omit_nl
:   if `TRUE` remove non-localized (NL) genes from the result.
Only applies at the gene or allele level.

sep
:   character defining both the input and output segment call 
delimiter.




Value
-------------------

A character vector containing allele, gene or family names.


References
-------------------

[http://imgt.org](http://imgt.org)



Examples
-------------------

```R
# Light chain examples
kappa_call <- c("Homsap IGKV1D-39*01 F,Homsap IGKV1-39*02 F,Homsap IGKV1-39*01",
"Homsap IGKJ5*01 F")

getAllele(kappa_call)

```


```
[1] "IGKV1-39*01" "IGKJ5*01"   

```


```R
getAllele(kappa_call, first=FALSE)

```


```
[1] "IGKV1-39*01,IGKV1-39*02" "IGKJ5*01"               

```


```R
getAllele(kappa_call, first=FALSE, strip_d=FALSE)

```


```
[1] "IGKV1D-39*01,IGKV1-39*02,IGKV1-39*01" "IGKJ5*01"                            

```


```R

getGene(kappa_call)

```


```
[1] "IGKV1-39" "IGKJ5"   

```


```R
getGene(kappa_call, first=FALSE)

```


```
[1] "IGKV1-39" "IGKJ5"   

```


```R
getGene(kappa_call, first=FALSE, strip_d=FALSE)

```


```
[1] "IGKV1D-39,IGKV1-39" "IGKJ5"             

```


```R

getFamily(kappa_call)

```


```
[1] "IGKV1" "IGKJ5"

```


```R
getFamily(kappa_call, first=FALSE)

```


```
[1] "IGKV1" "IGKJ5"

```


```R
getFamily(kappa_call, first=FALSE, collapse=FALSE)

```


```
[1] "IGKV1,IGKV1,IGKV1" "IGKJ5"            

```


```R
getFamily(kappa_call, first=FALSE, strip_d=FALSE)

```


```
[1] "IGKV1D,IGKV1" "IGKJ5"       

```


```R

# Heavy chain examples
heavy_call <- c("Homsap IGHV1-69*01 F,Homsap IGHV1-69D*01 F", 
"Homsap IGHD1-1*01 F", 
"Homsap IGHJ1*01 F")

getAllele(heavy_call, first=FALSE)

```


```
[1] "IGHV1-69*01" "IGHD1-1*01"  "IGHJ1*01"   

```


```R
getAllele(heavy_call, first=FALSE, strip_d=FALSE)

```


```
[1] "IGHV1-69*01,IGHV1-69D*01" "IGHD1-1*01"               "IGHJ1*01"                

```


```R

getGene(heavy_call, first=FALSE)

```


```
[1] "IGHV1-69" "IGHD1-1"  "IGHJ1"   

```


```R
getGene(heavy_call, first=FALSE, strip_d=FALSE)

```


```
[1] "IGHV1-69,IGHV1-69D" "IGHD1-1"            "IGHJ1"             

```


```R

# Filtering non-localized genes
nl_call <- c("IGHV3-NL1*01,IGHV3-30-3*01,IGHV3-30*01", 
"Homosap IGHV3-30*01 F,Homsap IGHV3-NL1*01 F",
"IGHV1-NL1*01")

getAllele(nl_call, first=FALSE, omit_nl=TRUE)

```


```
[1] "IGHV3-30-3*01,IGHV3-30*01" "IGHV3-30*01"               ""                         

```


```R
getGene(nl_call, first=FALSE, omit_nl=TRUE)

```


```
[1] "IGHV3-30-3,IGHV3-30" "IGHV3-30"            ""                   

```


```R
getFamily(nl_call, first=FALSE, omit_nl=TRUE)
```


```
[1] "IGHV3" "IGHV3" "IGHV1"

```



See also
-------------------

[countGenes](countGenes.md)



