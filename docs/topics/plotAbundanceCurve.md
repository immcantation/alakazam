**plotAbundanceCurve** - *Plots a clonal abundance distribution*


Usage
--------------------
```
plotAbundanceCurve(
data,
colors = NULL,
main_title = "Rank Abundance",
legend_title = NULL,
xlim = NULL,
ylim = NULL,
annotate = c("none", "depth"),
silent = FALSE,
...
)
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
`c(lower, upper)` x-axis limits. The lower x-axis 
value must be >=1.

ylim
:   numeric vector of two values specifying the 
`c(lower, upper)` y-axis limits. The limits on the 
abundance values are expressed as fractions of 1: use
c(0,1) to set the lower and upper limits to 0

\itemannotatestring defining whether to added values to the group labels 
of the legend. When `"none"` (default) is specified no
annotations are added. Specifying (`"depth"`) adds 
sequence counts to the labels.

\itemsilentif `TRUE` do not draw the plot and just return the ggplot2 
object; if `FALSE` draw the plot.

\item...additional arguments to pass to ggplot2::theme.












