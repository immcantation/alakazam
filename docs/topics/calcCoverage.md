





**calcCoverage** - *Calculate sample coverage*

Description
--------------------

`calcCoverage` calculates the sample coverage estimate, a measure of sample 
completeness, for varying orders using the method of Chao et al, 2015, falling back 
to the Chao1 method in the first order case.


Usage
--------------------
```
calcCoverage(x, r = 1)
```

Arguments
-------------------

x
:   numeric vector of abundance counts.

r
:   coverage order to calculate.



Value
-------------------

The sample coverage of the given order `r`.

References
-------------------


1. Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
Scand J Stat. 1984 11, 265270.
1. Chao A, et al. Unveiling the species-rank abundance distribution by 
generalizing the Good-Turing sample coverage theory. 
Ecology. 2015 96, 11891201.




Examples
-------------------

```R
# Load example data
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

# Calculate clone sizes
clones <- countClones(df, groups="SAMPLE")
# Calculate 1st order coverage for a single sample
calcCoverage(clones$SEQ_COUNT[clones$SAMPLE == "RL01"])
```


```
[1] 0.1608073

```



See also
-------------------

Used by [rarefyDiversity](rarefyDiversity.md).


