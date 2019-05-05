**sortGenes** - *Sort V(D)J genes*

Description
--------------------

`sortGenes` sorts a vector of V(D)J gene names by either lexicographic ordering 
or locus position.


Usage
--------------------
```
sortGenes(genes, method = c("name", "position"))
```

Arguments
-------------------

genes
:   vector of strings respresenting V(D)J gene names.

method
:   string defining the method to use for sorting genes. One of:

+  `"name"`:      sort in lexicographic order. Order is by 
family first, then gene, and then allele. 
+  `"position"`:  sort by position in the locus, as
determined by the final two numbers 
in the gene name. Non-localized genes 
are assigned to the highest positions.





Value
-------------------

A sorted character vector of gene names.



Examples
-------------------

```R
# Create a list of allele names
genes <- c("IGHV1-69D*01","IGHV1-69*01","IGHV4-38-2*01","IGHV1-69-2*01",
"IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
"IGHV1-2", "IGHV1-2*02", "IGHV1-69*02")

# Sort genes by name
sortGenes(genes)

```


```
 [1] "IGHV1-2"               "IGHV1-2*01,IGHV1-2*05" "IGHV1-2*02"           
 [4] "IGHV1-69*01"           "IGHV1-69D*01"          "IGHV1-69*02"          
 [7] "IGHV1-69-2*01"         "IGHV1-NL1*01"          "IGHV2-5*01"           
[10] "IGHV4-38-2*01"        

```


```R

# Sort genes by position in the locus
sortGenes(genes, method="pos")
```


```
 [1] "IGHV1-NL1*01"          "IGHV1-69-2*01"         "IGHV1-69*01"          
 [4] "IGHV1-69D*01"          "IGHV1-69*02"           "IGHV4-38-2*01"        
 [7] "IGHV2-5*01"            "IGHV1-2"               "IGHV1-2*01,IGHV1-2*05"
[10] "IGHV1-2*02"           

```



See also
-------------------

See `getAllele`, `getGene` and `getFamily` for parsing
gene names.



