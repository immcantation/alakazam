





**writeChangeoDb** - *Write a Change-O tab-delimited database file*

Description
--------------------

`writeChangeoDb` is a simple wrapper around [write_delim](http://www.rdocumentation.org/packages/readr/topics/write_delim) with defaults 
appropriate for writing a Change-O tab-delimited database file from a data.frame.


Usage
--------------------
```
writeChangeoDb(data, file)
```

Arguments
-------------------

data
:   data.frame of Change-O data.

file
:   output file name.





Examples
-------------------

```R
### Not run:
# Write a database
# writeChangeoDb(data, "changeo.tsv")
```



See also
-------------------

Wraps [write_delim](http://www.rdocumentation.org/packages/readr/topics/write_delim). See [readChangeoDb](readChangeoDb.md) for reading to Change-O files.



