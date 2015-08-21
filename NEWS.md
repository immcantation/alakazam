Version 0.2.1.beta:  August 21, 2015
-------------------------------------------------------------------------------

General:

+ Renamed `getDNADistMatrix()` and `getAADistMatrix()` to `getDNAMatrix` and 
  `getAAMatrix()`, respectively.
+ Added `getSeqMatrix()` which calculates a pairwise distance matrix for a set 
  of sequences.

Annotation:

+ Added support for unusual TCR gene names, such as 'TRGVA*01'.
+ Added removal of 'D' label (gene duplication) from gene names when parsed 
  with `getSegment()`, `getAllele()`, `getGene()` and `getFamily()`.  May be 
  disabled by providing the argument `strip_d=FALSE`.
  

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