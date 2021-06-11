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
[1] "/tmp/Rtmplu7Qeo/Clone50-temp-a656e656484"

```



See also
-------------------

This is just a wrapper for [tempfile](http://www.rdocumentation.org/packages/base/topics/tempfile) and 
[dir.create](http://www.rdocumentation.org/packages/base/topics/files2).






