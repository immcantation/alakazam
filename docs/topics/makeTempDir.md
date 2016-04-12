





**makeTempDir** - *Create a temporary folder*

Description
--------------------

`makeTempDir` creates a randomly named temporary folder in the 
system temp location.


Usage
--------------------
```
makeTempDir(prefix)
```

Arguments
-------------------

prefix
:   prefix name for the folder.



Value
-------------------

The path to the temporary folder.



Examples
-------------------

```R
makeTempDir("Clone50")
```


```
[1] "/tmp/RtmplOYXR4/Clone50-temp-371965744949"

```



See also
-------------------

This is just a wrapper for [tempfile](http://www.inside-r.org/r-doc/base/tempfile) and 
[dir.create](http://www.inside-r.org/r-doc/base/files2).



