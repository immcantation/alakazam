# Fastq

The `alakazam` package includes a set of functions to inspect the sequencing quality.

## Example data

Load example data:


```r
library(alakazam)
library(dplyr)
library(airr)

db <- read_rearrangement(system.file("extdata", "example_quality.tsv", package="alakazam"))
```

```
## Error: '' does not exist in current working directory ('/home/susanna/Documents/Work/Yale/projects/software_projects/alakazam/vignettes').
```

```r
fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")
```

## Load quality scores

This method allows to add the quality scores to the repertoire `data.frame` as strings.


```r
original_cols <- colnames(db)
db <- readFastqDb(db, fastq_file, style="both", quality_sequence=TRUE)
```

```
## Warning in read.FASTA(fl): failed to read sequences, returns NULL
```

```
## Error in attr(DNA, "QUAL") <- QUAL: attempt to set an attribute on NULL
```

```r
new_cols <- setdiff(colnames(db), original_cols)
db[,new_cols] %>% head()
```

```
## # A tibble: 6 x 0
```

The function `readFastq` takes as main inputs a repertoire `data.frame` (`db`) and 
a path to the corresponding `.fastq` file (`fastq_file`). The sequencing quality scores will
be merged into the `data.frame` by `sequence_id`. The newly added columns are:
. The other fields, contain the ASCII quality scores in the 
form of a vector, where values are comma separated, and `-` or `.` positions 
have value `" "` (blank).

After loading the quality scores with `readFastqDb`,  `getPositionQuality`
can be used to generate a `data.frame` of sequencing quality values 
per position.


```r
quality <- getPositionQuality(db, sequence_id="sequence_id", 
                              sequence="sequence_alignment",
                              quality_num="quality_alignment_num")
```

```
## Error in strsplit(data[[quality_num]][i], ","): non-character argument
```

```r
head(quality)
```

```
## Error in head(quality): object 'quality' not found
```


```r
min_pos <- min(quality$position)
```

```
## Error in eval(expr, envir, enclos): object 'quality' not found
```

```r
max_pos <- max(quality$position)
```

```
## Error in eval(expr, envir, enclos): object 'quality' not found
```

```r
ggplot(quality, aes(x=position,
                    y=quality_alignment_num,
                    color=nt)) +
  geom_point() +
  coord_cartesian(xlim=c(110,120)) +
  xlab("IMGT position") +
  ylab("Sequencing quality") +
  scale_fill_gradient(low = "light blue",  high = "dark red") +
  scale_x_continuous(breaks=c(min_pos:max_pos)) +
  alakazam::baseTheme()
```

```
## Error in ggplot(quality, aes(x = position, y = quality_alignment_num, : object 'quality' not found
```

You can add use the quality `data.frame` to complement analysis performed
with other tools from the Immcantation framework. For example, you could inspect
the sequencing quality of novel polymorphisms identified with `tigger`, or
the sequencing quality in mutated/unmutated regions.

## Mask low quality positions

Use `maskPositionsByQuality` to mask low quality positions. Positions with
a sequencing quality < `min_quality` will be replaced with an 'N'. A message
will show the number of sequences in `db` that had at least one position
masked.


```r
db <- maskPositionsByQuality(db, min_quality=70,
                             sequence="sequence_alignment",
                             quality="quality_alignment_num")
```

```
## Error in strsplit(db_row[[quality_num]], ","): non-character argument
```

