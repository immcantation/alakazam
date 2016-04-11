





**writeChangeoDb** - *Write a Change-O tab-delimited database file*

Description
--------------------

`writeChangeoDb` is a simple wrapper around [write.table](http://www.inside-r.org/r-doc/utils/write.table) with defaults 
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
# writeChangeoDb(data, "changeo_output.tab")
```



See also
-------------------

Wraps [write.table](http://www.inside-r.org/r-doc/utils/write.table). See [readChangeoDb](readChangeoDb.md) for reading to Change-O files.



