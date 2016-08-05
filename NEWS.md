Version 0.2.4.999:  August 5, 2016
-------------------------------------------------------------------------------

General:

+ Added stop and gap characters to `getAAMatrix()` return matrix.


Version 0.2.4:  July 20, 2016
-------------------------------------------------------------------------------

General:

+ Added Rcpp and data.table dependencies.
+ Modified `readChangeoDb()` to wrap `data.table::fread()` instead of 
  `utils::read.table()` if the input file is not compressed.
+ Ported `testSeqEqual()`, `getSeqDistance()` and `getSeqMatrix()` to C++ to 
  improve performance of `collapseDuplicates()` and other dependent functions.
+ Renamed `testSeqEqual()`, `getSeqDistance()` and `getSeqMatrix()` to 
  `seqEqual()`, `seqDist()` and `pairwiseDist()`, respectively.
+ Added `pairwiseEqual()` which creates a logical sequence distance matrix;
  TRUE if sequences are identical, FALSE if not, excluding Ns and gaps.
+ Added translation of ambiguous and gap characters to `X` in 
  `translateDNA()`.
+ Fixed bug in `collapseDuplicates()` wherein the input data type sanity check
  would cause the vignette to fail to build under R 3.3.
+ Replaced the `ExampleDb.gz` file with a larger, more clonal, `ExampleDb` 
  data object.
+ Replaced `ExampleTrees` with a larger set of trees.
+ Renamed `multiggplot()` to `gridPlot()`.

Amino Acid Analysis:

+ Set default to `normalize=FALSE` for charge calculations to be more consistent
  with previously published repertoire sequencing results.
  
Diversity Analysis:

+ Added a `progress` argument to `rarefyDiversity()` and `testDiversity()` to
  enable the (previously default) progress bar.
+ Fixed a bug in `estimateAbundance()` were the function would fail if there 
  was only a single input sequence per group.
+ Changed column names in `data` and `summary` slots of `DiversityTest` to 
  uppercase for consistency with other tools.
+ Added dispatching of `plot` to `plotDiversityCurve` for `DiversityCurve`
  objects.
  
Gene Usage:

+ Added `sortGenes()` function to sort V(D)J genes by name or locus position.
+ Added `clone` argument to `countGenes()` to allow restriction of gene 
  abundance to one gene per clone.

Topology Analysis:

+ Added a set of functions for lineage tree topology analysis.
+ Added a vignette showing basic tree topology analysis.


Version 0.2.3:  February 22, 2016
-------------------------------------------------------------------------------

General:

+ Fixed a bug wherein the package would not build on R < 3.2.0 due to changes
  in `base::nchar()`.
+ Changed R dependency to R >= 3.1.2.


Version 0.2.2:  January 29, 2016
-------------------------------------------------------------------------------

General:

+ Updated license from CC BY-NC-SA 3.0 to CC BY-NC-SA 4.0.
+ Internal changes to conform to CRAN policies.

Amino Acid Analysis:

+ Fixed bug where arguments for the `aliphatic()` function were not being
  passed through the ellipsis argument of `aminoAcidProperties()`.
+ Improved amino acid analysis vignette.
+ Added check for correctness of amino acids sequences to `aminoAcidProperties()`.
+ Renamed `AA_TRANS` to `ABBREV_AA`.

Diversity:

+ Added evenness and bootstrap standard deviation to `rarefyDiversity()` 
  output.

Lineage:

+ Added `ExampleTrees` data with example output from `buildPhylipLineage()`.


Version 0.2.1:  December 18, 2015
-------------------------------------------------------------------------------

General:

+ Removed plyr dependency.
+ Added dplyr, lazyeval and stringi dependencies.
+ Added strict requirement for igraph version >= 1.0.0.
+ Renamed `getDNADistMatrix()` and `getAADistMatrix()` to `getDNAMatrix` and 
  `getAAMatrix()`, respectively.
+ Added `getSeqMatrix()` which calculates a pairwise distance matrix for a set 
  of sequences.
+ Modified default plot sizing to be more appropriate for export to PDF 
  figures with 7-8 inch width.
+ Added `multiggplot()` function for performing multiple panel plots.

Amino Acid Analysis:

+ Migrated amino acid property analysis from Change-O CTL to alakazam. 
  Includes the new functions `gravy()`, `bulk()`, `aliphatic()`, `polar()`, 
  `charge()`, `countPatterns()` and `aminoAcidProperties()`.

Annotation:

+ Added support for unusual TCR gene names, such as 'TRGVA*01'.
+ Added removal of 'D' label (gene duplication) from gene names when parsed 
  with `getSegment()`, `getAllele()`, `getGene()` and `getFamily()`.  May be 
  disabled by providing the argument `strip_d=FALSE`.
+ Added `countGenes()` to tabulate V(D)J allele, gene and family usage.

Diversity:

+ Added several functions related to analysis of clone size distributions, 
  including `countClones()`, `estimateAbundance()` and `plotAbundance()`.
+ Renamed `resampleDiversity()` to `rarefyDiversity()` and changed many of
  the internals. Bootstrapping is now performed on an inferred complete
  relative abundance distribution.
+ Added support for inclusion of copy number in clone size determination
  within `rarefyDiversity()` and `testDiversity()`.
+ Diversity scores and confiderence intervals within `rarefyDiversity()`
  and `testDiversity()` are now calculated using the mean and standard 
  deviation of the bootstrap realizations, rather than the median and
  upper/lower quantiles.
+ Added ability to add counts to the legend in `plotDiversityCurve()`.


Version 0.2.0:  June 15, 2015
-------------------------------------------------------------------------------

Initial public release.

General:

+ Added citations for the `citation("alakazam")` command.


Version 0.2.0.beta-2015-05-30:  May 30, 2015
-------------------------------------------------------------------------------

Lineage:

+ Added more error checking to `buildPhylipLineage()`.


Version 0.2.0.beta-2015-05-26:  May 26, 2015
-------------------------------------------------------------------------------

Lineage:

+ Fixed issue where `buildPhylipLineage()` would hang on R 3.2 due to R change 
  request PR#15508.


Version 0.2.0.beta-2015-05-05:  May 05, 2015
-------------------------------------------------------------------------------

Prerelease for review.
