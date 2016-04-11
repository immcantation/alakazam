





**countGenes** - *Tabulates V(D)J allele, gene or family usage.*

Description
--------------------

Determines the count and relative abundance of V(D)J alleles, genes or families within
groups.

Usage
--------------------

```
countGenes(data, gene, groups = NULL, copy = NULL, mode = c("gene",
"allele", "family"))
```

Arguments
-------------------

data
:   data.frame with Change-O style columns containing clonal assignments.

gene
:   column containing allele assignments. Only the first allele in the
column will be considered.

groups
:   columns containing grouping variables. If `NULL` do not group.

copy
:   name of the `data` column containing copy numbers for each 
sequence. If this value is specified, then total copy abundance
is determined by the sum of copy numbers within each gene.

mode
:   one of `c("gene", "family", "allele")` defining
the degree of specificity regarding allele calls. Determines whether 
to return counts for genes, families or alleles.



Value
-------------------

A data.frame summarizing family, gene or allele counts and frequencies with
columns:

+  `GENE`:        name of the family, gene or allele
+  `SEQ_COUNT`:   total number of sequences for the gene.
+  `SEQ_FREQ`:    frequency of the gene as a fraction of the total
number of sequences within each grouping.
+  `COPY_COUNT`:  sum of the copy counts in the `copy` column.
for each gene. Only present if the `copy` 
argument is specified.
+  `COPY_FREQ`:   frequency of the gene as a fraction of the total
copy number within each group. Only present if 
the `copy` argument is specified.

Additional columns defined by the `groups` argument will also be present.



Examples
-------------------

```R
# Load example data
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

# Without copy numbers
genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="family")
genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="gene")
genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="allele")

# With copy numbers and multiple groups
genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
copy="DUPCOUNT", mode="family")
genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
copy="DUPCOUNT", mode="gene")
genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
copy="DUPCOUNT", mode="allele")
```




