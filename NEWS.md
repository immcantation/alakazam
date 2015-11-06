Version 0.2.1.beta:  October 24, 2015
-------------------------------------------------------------------------------

General:

+ Removed plyr dependency. Added dplyr, iNEXT, lazyeval and stringi 
  dependencies.
+ Renamed `getDNADistMatrix()` and `getAADistMatrix()` to `getDNAMatrix` and 
  `getAAMatrix()`, respectively.
+ Added `getSeqMatrix()` which calculates a pairwise distance matrix for a set 
  of sequences.

Annotation:

+ Added support for unusual TCR gene names, such as 'TRGVA*01'.
+ Added removal of 'D' label (gene duplication) from gene names when parsed 
  with `getSegment()`, `getAllele()`, `getGene()` and `getFamily()`.  May be 
  disabled by providing the argument `strip_d=FALSE`.
+ Added `countGenes()` to tabulate V(D)J allele, gene and family usage.

Clonality:

+ Added several functions related to analysis of clone size distributions, 
  including `countClones()`, `estimateAbundance()` and `plotAbundance`.
  
Diversity:

+ Renamed `resampleDiversity()` to `rarefyDiversity()` and changed many of
  the internals. Bootstrapping is now performed on an inferred complete
  relative abundance distribution, using methods from the iNEXT package.
+ Added support for inclusion of copy number in clone size determination
  within `rarefyDiversity()` and `testDiversity()`.
+ Added coverage based rarefaction and no rarefaction to the options for 
  resampling strategies in `rarefyDiversity()`.
+ Diversity scores and confiderence intervals within `rarefyDiversity()`
  and `testDiversity()` are now calculated using the mean and standard 
  deviation of the bootstrap realizations, rather than the median and
  upper/lower quantiles.
+ Added ability to add counts and coverage values to the legend in
  `potDiversityCurve()`.


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