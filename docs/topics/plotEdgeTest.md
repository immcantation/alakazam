**plotEdgeTest** - *Plot the results of an edge permutation test*

Description
--------------------

`plotEdgeTest` plots the results of an edge permutation test performed with 
`testEdges` as either a histogram or cumulative distribution function.


Usage
--------------------
```
plotEdgeTest(data, color = "black", main_title = "Edge Test",
style = c("histogram", "cdf"), silent = FALSE, ...)
```

Arguments
-------------------

data
:   [EdgeTest](EdgeTest-class.md) object returned by [testEdges](testEdges.md).

color
:   color of the histogram or lines.

main_title
:   string specifying the plot title.

style
:   type of plot to draw. One of:

+  `"histogram"`:  histogram of the edge count 
distribution with a red dotted line
denoting the observed value.
+  `"cdf"`:        cumulative distribution function 
of edge counts with a red dotted 
line denoting the observed value and
a blue dotted line indicating the 
p-value.


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
# Define example tree set
graphs <- ExampleTrees[1-10]

# Perform edge test on isotypes
x <- testEdges(graphs, "ISOTYPE", nperm=10)

```

**Error in ecdf(d)**: 'x' must have 1 or more non-missing values
```R

# Plot
plotEdgeTest(x, color="steelblue", style="hist")

```

**Error in rename(data@tests, Parent = "PARENT", Child = "CHILD")**: object 'x' not found
```R
plotEdgeTest(x, style="cdf")
```

**Error in rename(data@tests, Parent = "PARENT", Child = "CHILD")**: object 'x' not found

See also
-------------------

See [testEdges](testEdges.md) for performing the test.



