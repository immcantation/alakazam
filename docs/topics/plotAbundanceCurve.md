**plotAbundanceCurve** - *Plots a clonal abundance distribution*

Description
--------------------

`plotAbundanceCurve` plots the results from estimating the complete clonal 
relative abundance distribution. The distribution is plotted as a log rank abundance 
distribution.


Usage
--------------------
```
plotAbundanceCurve(data, colors = NULL, main_title = "Rank Abundance",
legend_title = NULL, xlim = NULL, ylim = NULL,
annotate = c("none", "depth"), silent = FALSE, ...)
```

Arguments
-------------------

data
:   [AbundanceCurve](AbundanceCurve-class.md) object returned by [estimateAbundance](estimateAbundance.md).

colors
:   named character vector whose names are values in the 
`group` column of `data` and whose values are 
colors to assign to those group names.

main_title
:   string specifying the plot title.

legend_title
:   string specifying the legend title.

xlim
:   numeric vector of two values specifying the 
`c(lower, upper)` x-axis limits.

ylim
:   numeric vector of two values specifying the 
`c(lower, upper)` y-axis limits.

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
# Estimate abundance by sample and plot
abund <- estimateAbundance(ExampleDb, group="SAMPLE", nboot=100)
plotAbundanceCurve(abund, legend_title="Sample")
```

**Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, **: polygon edge not found![3](plotAbundanceCurve-3.png)


See also
-------------------

See [AbundanceCurve](AbundanceCurve-class.md) for the input object and [estimateAbundance](estimateAbundance.md) for 
generating the input abundance distribution.
Plotting is performed with [ggplot](http://www.rdocumentation.org/packages/ggplot2/topics/ggplot).



