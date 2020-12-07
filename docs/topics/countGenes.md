**countGenes** - *Tabulates V(D)J allele, gene or family usage.*

Description
--------------------

Determines the count and relative abundance of V(D)J alleles, genes or families within
groups.


Usage
--------------------
```
countGenes(
data,
gene,
groups = NULL,
copy = NULL,
clone = NULL,
fill = FALSE,
mode = c("gene", "allele", "family", "asis"),
remove_na = TRUE
)
```

Arguments
-------------------

data
:   data.frame with AIRR-format or Change-O style columns.

gene
:   column containing allele assignments. Only the first allele in the
column will be considered when `mode` is "gene", "family" or 
"allele". The value will be used as it is with `mode="asis"`.

groups
:   columns containing grouping variables. If `NULL` do not group.

copy
:   name of the `data` column containing copy numbers for each 
sequence. If this value is specified, then total copy abundance
is determined by the sum of copy numbers within each gene.
This argument is ignored if `clone` is specified.

clone
:   name of the `data` column containing clone identifiers for each 
sequence. If this value is specified, then one gene will be considered 
for each clone. Note, this is accomplished by using the most 
common gene within each `clone` identifier. As such,
ambiguous alleles within a clone will not be accurately represented.

fill
:   logical of `c(TRUE, FALSE)` specifying when if groups (when specified)
lacking a particular gene should be counted as 0 if TRUE or not (omitted)

mode
:   one of `c("gene", "family", "allele", "asis")` defining
the degree of specificity regarding allele calls. Determines whether 
to return counts for genes (calling `getGene`), 
families (calling `getFamily`), alleles (calling 
`getAllele`) or using the value as it is in the column
`gene`, without any processing.

remove_na
:   removes rows with `NA` values in the gene column if `TRUE` and issues a warning. 
Otherwise, keeps those rows and considers `NA` as a gene in the final counts 
and relative abundances.




Value
-------------------

A data.frame summarizing family, gene or allele counts and frequencies 
with columns:

+  `gene`:         name of the family, gene or allele
+  `seq_count`:    total number of sequences for the gene.
+  `seq_freq`:     frequency of the gene as a fraction of the total
number of sequences within each grouping.
+  `copy_count`:   sum of the copy counts in the `copy` column.
for each gene. Only present if the `copy` 
argument is specified.
+  `copy_freq`:    frequency of the gene as a fraction of the total
copy number within each group. Only present if 
the `copy` argument is specified.
+  `clone_count`:  total number of clones for the gene.
+  `clone_freq`:   frequency of the gene as a fraction of the total
number of clones within each grouping.

Additional columns defined by the `groups` argument will also be present.



Examples
-------------------

```R
# Without copy numbers
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="family")
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="gene")
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="allele")

# With copy numbers and multiple groups
genes <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
copy="duplicate_count", mode="family")

# Count by clone
genes <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
clone="clone_id", mode="family")

# Count absent genes 
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", 
mode="allele", fill=TRUE)
```








