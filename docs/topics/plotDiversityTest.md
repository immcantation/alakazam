





**plotDiversityTest** - *Plot the results of TestDiversity*

Description
--------------------

`plotDiversityTest` plots a `DiversityTest` object as the mean
with a line range indicating plus/minus one standard deviation.


Usage
--------------------
```
plotDiversityTest(data, colors = NULL, main_title = "Diversity",
legend_title = "Group", log_d = FALSE, annotate = c("none", "depth"),
silent = FALSE, ...)
```

Arguments
-------------------

data
:   [DiversityTest](DiversityTest-class.md) object returned by 
[testDiversity](testDiversity.md).

colors
:   named character vector whose names are values in the 
`group` column of the `data` slot of `data`,
and whose values are colors to assign to those group names.

main_title
:   string specifying the plot title.

legend_title
:   string specifying the legend title.

log_d
:   if `TRUE` then plot the diversity scores <code class = 'eq'>D</code> 
on a log scale; if `FALSE` plot on a linear scale.

annotate
:   string defining whether to added values to the group labels 
of the legend. When `"none"` (default) is specified no
annotations are added. Specifying (`"depth"`) adds 
sequence counts to the labels.

silent
:   if `TRUE` do not draw the plot and just return the ggplot2 
object; if `FALSE` draw the plot.

...
:   additional arguments to pass to ggplot2::theme.




Value
-------------------

A `ggplot` object defining the plot.



Examples
-------------------

```R
# All groups pass default minimum sampling threshold of 10 sequences
div <- testDiversity(ExampleDb, "SAMPLE", q=0, nboot=100)
plotDiversityTest(div, legend_title="Sample")
```

![2](plotDiversityTest-2.png)


See also
-------------------

See [testDiversity](testDiversity.md) for generating [DiversityTest](DiversityTest-class.md)
objects for input. Plotting is performed with [ggplot](http://www.rdocumentation.org/packages/ggplot2/topics/ggplot).



