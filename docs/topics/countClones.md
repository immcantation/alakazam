**countClones** - *Tabulates clones sizes*

Description
--------------------

`countClones` determines the number of sequences and total copy number of 
clonal groups.


Usage
--------------------
```
countClones(data, groups = NULL, copy = NULL, clone = "clone_id")
```

Arguments
-------------------

data
:   data.frame with Change-O style columns containing clonal assignments.

groups
:   character vector defining `data` columns containing grouping 
variables. If `groups=NULL`, then do not group data.

copy
:   name of the `data` column containing copy numbers for each 
sequence. If this value is specified, then total copy abundance
is determined by the sum of copy numbers within each clonal group.

clone
:   name of the `data` column containing clone identifiers.




Value
-------------------

A data.frame summarizing clone counts and frequencies with columns:

+  `CLONE`:    	 clone identifier.
+  `SEQ_COUNT`:   total number of sequences for the clone.
+  `SEQ_FREQ`:    frequency of the clone as a fraction of the total
number of sequences within each group.
+  `COPY_COUNT`:  sum of the copy counts in the `copy` column.
Only present if the `copy` argument is 
specified.
+  `COPY_FREQ`:   frequency of the clone as a fraction of the total
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








