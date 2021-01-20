**countClones** - *Tabulates clones sizes*

Description
--------------------

`countClones` determines the number of sequences and total copy number of 
clonal groups.


Usage
--------------------
```
countClones(
data,
groups = NULL,
copy = NULL,
clone = "clone_id",
remove_na = TRUE
)
```

Arguments
-------------------

data
:   data.frame with columns containing clonal assignments.

groups
:   character vector defining `data` columns containing grouping 
variables. If `groups=NULL`, then do not group data.

copy
:   name of the `data` column containing copy numbers for each 
sequence. If this value is specified, then total copy abundance
is determined by the sum of copy numbers within each clonal group.

clone
:   name of the `data` column containing clone identifiers.

remove_na
:   removes rows with `NA` values in the clone column if `TRUE` and issues a warning. 
Otherwise, keeps those rows and considers `NA` as a clone in the final counts 
and relative abundances.




Value
-------------------

A data.frame summarizing clone counts and frequencies with columns:

+  `clone_id`:    clone identifier. This is the default column
name, specified with `clone='clone_id'`.
If the function call uses Change-O 
formatted data and `clone='CLONE'`, this
column will have name `CLONE`.
+  `seq_count`:   total number of sequences for the clone.
+  `seq_freq`:    frequency of the clone as a fraction of the total
number of sequences within each group.
+  `copy_count`:  sum of the copy counts in the `copy` column.
Only present if the `copy` argument is 
specified.
+  `copy_freq`:   frequency of the clone as a fraction of the total
copy number within each group. Only present if 
the `copy` argument is specified.

Also includes additional columns specified in the `groups` argument.



Examples
-------------------

```R
# Without copy numbers
clones <- countClones(ExampleDb, groups="sample_id")

# With copy numbers and multiple groups
clones <- countClones(ExampleDb, groups=c("sample_id", "c_call"), copy="duplicate_count")
```








