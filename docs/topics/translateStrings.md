





**translateStrings** - *Translate a vector of strings*

Description
--------------------

`translateStrings` modifies a character vector by substituting one or more 
strings with a replacement string.

Usage
--------------------

```
translateStrings(strings, translation)
```

Arguments
-------------------

strings
:   vector of character strings to modify.

translation
:   named character vector or a list of character vectors specifying 
the strings to replace (values) and their replacements (names).



Value
-------------------

A modified `strings` vector.

Details
-------------------

Does not perform partial replacements. Each translation value must match a complete 
`strings` value or it will not be replaced.  Values that do not have a replacement
named in the `translation` parameter will not be modified.

Replacement is accomplished using [gsub](http://www.inside-r.org/r-doc/base/grep).



Examples
-------------------

```R
# Using a vector translation
strings <- LETTERS[1:5]
translation <- c("POSITION1"="A", "POSITION5"="E")
translateStrings(strings, translation)

```


```
[1] "POSITION1" "B"         "C"         "D"         "POSITION5"

```


```R

# Using a list translation
strings <- LETTERS[1:5]
translation <- list("1-3"=c("A","B","C"), "4-5"=c("D","E"))
translateStrings(strings, translation)
```


```
[1] "1-3" "1-3" "1-3" "4-5" "4-5"

```



See also
-------------------

See [gsub](http://www.inside-r.org/r-doc/base/grep) for single value replacement in the base package.



