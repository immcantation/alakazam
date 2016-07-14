Diversity analysis
====================


The clonal diversity of the repertoire can be analyzed using the general form
of the diversity index, as proposed by Hill in:

    Hill, M. Diversity and evenness: a unifying notation and its consequences. 
        Ecology 54, 427-432 (1973).

Coupled with resampling strategies to correct for variations in sequencing 
depth, as well as inferrence of complete clonal abundance distributions as 
described in:

    Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
        A framework for sampling and estimation in species diversity studies. 
        Ecol Monogr. 2014 84:45-67.
    Chao A, et al. Unveiling the species-rank abundance distribution by 
        generalizing the Good-Turing sample coverage theory. 
        Ecology. 2015 96, 11891201.

This package provides methods for the inferrence of a complete clonal 
abundance distribution, using the `estimateAbundance` function, along with 
two approaches to assess diversity of these distributions: 

1. Generation of a smooth diversity (D) curve over a range of diversity orders (q) 
using `rarefyDiversity`.
2. A significance test of the diversity (D) at a fixed diversity order (q) using 
`testDiversity`.


## Load Change-O data

A small example Change-O tab-delimited database file is included in the 
`alakazam` package. Diversity calculation requires the `CLONE` field 
(column) to be present in the Change-O file, as well as an additional grouping 
column. In this example we will use the grouping columns `SAMPLE` and `ISOTYPE`.


```r
library(alakazam)

# Load Change-O file
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
db <- readChangeoDb(file)
```

## Generate a clonal abundance curve

A simple table of the observed clonal abundance counts and frequencies may be
generated using the `countClones` function either without copy numbers, where
the size of each clone is determined by the number of sequence members:


```r
# Partitions the data based on the SAMPLE column
clones <- countClones(db, groups="SAMPLE")
head(clones, 5)
```

```
## Source: local data frame [5 x 4]
## Groups: SAMPLE [1]
## 
##   SAMPLE CLONE SEQ_COUNT SEQ_FREQ
##    <chr> <chr>     <int>    <dbl>
## 1    +7d  3128       100    0.100
## 2    +7d  3100        50    0.050
## 3    +7d  3141        44    0.044
## 4    +7d  3177        30    0.030
## 5    +7d  3170        28    0.028
```

You may also specify a column containing the abundance count of each sequence 
(usually copy numbers), that will including weighting of each clone size by the 
corresponding abundance count. Furthermore, multiple grouping columns may be
specified such that `SEQ_FREQ` (unwieghted clone size as a fraction
of total sequences in the group) and `COPY_FREQ` (weighted faction) are 
normalized to within multiple group data partitions.


```r
# Partitions the data based on both the SAMPLE and ISOTYPE columns
# Weights the clone sizes by the DUPCOUNT column
clones <- countClones(db, groups=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")
head(clones, 5)
```

```
## Source: local data frame [5 x 7]
## Groups: SAMPLE, ISOTYPE [2]
## 
##   SAMPLE ISOTYPE CLONE SEQ_COUNT COPY_COUNT   SEQ_FREQ  COPY_FREQ
##    <chr>   <chr> <chr>     <int>      <int>      <dbl>      <dbl>
## 1    +7d     IgA  3128        88        651 0.33082707 0.49732620
## 2    +7d     IgG  3100        49        279 0.09280303 0.17296962
## 3    +7d     IgA  3141        44        240 0.16541353 0.18334607
## 4    +7d     IgG  3192        19        141 0.03598485 0.08741476
## 5    +7d     IgG  3177        29        130 0.05492424 0.08059516
```

While `countClones` will report observed abundances, it will not correct the
distribution nor provide confidence intervals. A complete clonal abundance 
distribution may be inferred using the `estimateAbundance` function with
confidence intervals derived via bootstrapping.  This output may be visualized
using the `plotAbundance` function.


```r
# Partitions the data on the SAMPLE column
# Calculates a 95% confidence interval via 200 bootstrap realizations
clones <- estimateAbundance(db, group="SAMPLE", ci=0.95, nboot=200)
```

```
## Error in estimateAbundance(db, group = "SAMPLE", ci = 0.95, nboot = 200): object 'progress' not found
```

```r
head(clones, 5)
```

```
## Source: local data frame [5 x 7]
## Groups: SAMPLE, ISOTYPE [2]
## 
##   SAMPLE ISOTYPE CLONE SEQ_COUNT COPY_COUNT   SEQ_FREQ  COPY_FREQ
##    <chr>   <chr> <chr>     <int>      <int>      <dbl>      <dbl>
## 1    +7d     IgA  3128        88        651 0.33082707 0.49732620
## 2    +7d     IgG  3100        49        279 0.09280303 0.17296962
## 3    +7d     IgA  3141        44        240 0.16541353 0.18334607
## 4    +7d     IgG  3192        19        141 0.03598485 0.08741476
## 5    +7d     IgG  3177        29        130 0.05492424 0.08059516
```

```r
# Plots a rank abundance curve of the relative clonal abundances
p1 <- plotAbundance(clones, legend_title="Sample")
```

```
## Error in eval(expr, envir, enclos): object 'LOWER' not found
```

![plot of chunk Diversity-Vignette-5](figure/Diversity-Vignette-5-1.png)


## Generate a diversity curve

The function `rarefyDiversity` performs uniform resampling of the input 
sequences and recalculates the clone size distribution, and diversity, with each 
resampling realization. Diversity (D) is calculated over a range of diversity 
orders (q) to generate a smooth curve.


```r
# Compare diversity curve across values in the "SAMPLE" column
# q ranges from 0 (min_q=0) to 32 (max_q=32) in 0.05 incriments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 2000 resampling realizations are performed (nboot=200)
sample_div <- rarefyDiversity(db, "SAMPLE", min_q=0, max_q=32, step_q=0.05, 
                                 ci=0.95, nboot=200)

# Compare diversity curve across values in the "ISOTYPE" column
# Analyse is restricted to ISOTYPE values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
isotype_div <- rarefyDiversity(db, "ISOTYPE", min_n=30, min_q=0, max_q=32, 
                                  step_q=0.05, ci=0.95, nboot=200)
```


```r
# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")
p2 <- plotDiversityCurve(sample_div, main_title=sample_main, 
                         legend_title="Sample", log_q=TRUE, log_d=TRUE)
```

![plot of chunk Diversity-Vignette-7](figure/Diversity-Vignette-7-1.png)

```r
# Plot isotype diversity using default set of Ig isotype colors
isotype_main <- paste0("Isotype diversity (n=", isotype_div@n, ")")
p3 <- plotDiversityCurve(isotype_div, colors=IG_COLORS, main_title=isotype_main, 
                         legend_title="Isotype", log_q=TRUE, log_d=TRUE)
```

![plot of chunk Diversity-Vignette-7](figure/Diversity-Vignette-7-2.png)

## Test diversity at a fixed diversity order

The function `testDiversity` performs resampling and diversity calculation in 
the same manner as `rarefyDiversity`, but only for a single diversity order. 
Significance testing across groups is performed using the delta of the bootstrap
distributions between groups.


```r
# Test diversity at q=0 (species richness) across values in the "SAMPLE" column
# 2000 bootstrap realizations are performed (nboot=200)
sample_test <- testDiversity(db, 0, "SAMPLE", nboot=200)
```

```r
sample_test
```

```
## An object of class "DiversityTest"
## Slot "tests":
##         test pvalue delta_mean delta_sd
## 1 -1h != +7d      0    480.425 16.77083
## 
## Slot "summary":
##     group    mean       sd
## -1h   -1h 818.965 12.02776
## +7d   +7d 338.540 12.77161
## 
## Slot "groups":
## [1] "-1h" "+7d"
## 
## Slot "q":
## [1] 0
## 
## Slot "n":
##  -1h  +7d 
## 1000 1000 
## 
## Slot "nboot":
## [1] 200
```

```r
# Test diversity across values in the "ISOTYPE" column
# Analyse is restricted to ISOTYPE values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
isotype_test <- testDiversity(db, 2, "ISOTYPE", min_n=30, nboot=200)
```

```r
isotype_test
```

```
## An object of class "DiversityTest"
## Slot "tests":
##         test pvalue delta_mean  delta_sd
## 1 IgA != IgD      0  191.98212 10.895729
## 2 IgA != IgG      0   26.82262  4.349422
## 3 IgA != IgM      0  228.85129  6.535732
## 4 IgD != IgG      0  165.15950 11.621868
## 5 IgD != IgM      0   36.86917 12.121182
## 6 IgG != IgM      0  202.02868  7.257912
## 
## Slot "summary":
##     group      mean        sd
## IgA   IgA  12.79637  1.990933
## IgD   IgD 204.77848 10.583917
## IgG   IgG  39.61898  4.200877
## IgM   IgM 241.64766  6.493150
## 
## Slot "groups":
## [1] "IgA" "IgD" "IgG" "IgM"
## 
## Slot "q":
## [1] 2
## 
## Slot "n":
## IgA IgD IgG IgM 
## 259 259 259 259 
## 
## Slot "nboot":
## [1] 200
```
