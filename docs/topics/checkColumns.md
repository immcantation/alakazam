**checkColumns** - *Check data.frame for valid columns and issue message if invalid*

Description
--------------------

Check data.frame for valid columns and issue message if invalid


Usage
--------------------
```
checkColumns(data, columns, logic = c("all", "any"))
```

Arguments
-------------------

data
:   data.frame to check.

columns
:   vector of column names to check.

logic
:   one of `"all"` or `"any"` controlling whether all,
or at least one, of the columns must be valid, respectively.




Value
-------------------

`TRUE` if columns are valid and a string message if not.



Examples
-------------------

```R
df <- data.frame(A=1:3, B=4:6, C=rep(NA, 3))
checkColumns(df, c("A", "B"), logic="all")

```


```
[1] TRUE

```


```R
checkColumns(df, c("A", "B"), logic="any")

```


```
[1] TRUE

```


```R
checkColumns(df, c("A", "C"), logic="all")

```


```
[1] "The column C contains no data"

```


```R
checkColumns(df, c("A", "C"), logic="any")

```


```
[1] TRUE

```


```R
checkColumns(df, c("A", "D"), logic="all")

```


```
[1] "The column D was not found"

```


```R
checkColumns(df, c("A", "D"), logic="any")
```


```
[1] TRUE

```








