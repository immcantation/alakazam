





**alakazam** - *The alakazam package*

Description
--------------------

`alakazam` in a member of the Change-O suite of tools and serves five main 
purposes:

+ Providing core functionality for other R packages in the Change-O suite. This
includes common tasks such as file I/O, basic DNA sequence manipulation, and
interacting with V(D)J segment and gene annotations.
+ Providing an R interface for interacting with the output of the pRESTO 
tool suite.
+ Performing lineage reconstruction on clonal populations of immunoglobulin 
(Ig) sequences. 
+ Performing clonal abundance and diversity analysis on lymphocyte repertoires.
+ Performing physicochemical property analyses of lymphocyte receptor sequences.

For additional details regarding the use of the `alakazam` package see the 
vignettes:
`browseVignettes("alakazam")`



File I/O
-------------------



+ [readChangeoDb](readChangeoDb.md):        Input Change-O style files.
+ [writeChangeoDb](writeChangeoDb.md):       Output Change-O style files.


Sequence cleaning
-------------------



+ [maskSeqEnds](maskSeqEnds.md):          Mask ragged ends.
+ [maskSeqGaps](maskSeqGaps.md):          Mask gap characters.
+ [collapseDuplicates](collapseDuplicates.md):   Remove duplicate sequences.


Lineage reconstruction
-------------------



+ [makeChangeoClone](makeChangeoClone.md):     Clean sequences for lineage reconstruction.
+ [buildPhylipLineage](buildPhylipLineage.md):   Perform lineage reconstruction of Ig sequences.


Diversity analysis
-------------------



+ [countClones](countClones.md):          Calculate clonal abundance.
+ [estimateAbundance](estimateAbundance.md):    Infer complete clonal abundance distribution with
confidence intervals.
+ [rarefyDiversity](rarefyDiversity.md):      Generate clonal diversity curves.
+ [testDiversity](testDiversity.md):        Test significance of clonal diversity scores.
+ [plotAbundance](plotAbundance.md):        Plot clone size distribution as a rank-abundance 
curve.
+ [plotDiversityCurve](plotDiversityCurve.md):   Plot clonal diversity curves.


Ig and TCR sequence annotation
-------------------



+ [countGenes](countGenes.md):           Calculate Ig and TCR allele, gene and family usage.
+ [extractVRegion](extractVRegion.md):       Extract CDRs and FWRs sub-sequences.
+ [getAllele](getSegment.md):            Get V(D)J allele names.
+ [getGene](getSegment.md):              Get V(D)J gene names.
+ [getFamily](getSegment.md):            Get V(D)J family names.


Sequence distance calculation
-------------------



+ [getSeqDistance](getSeqDistance.md):       Calculate Hamming distance between sequences.
+ [getSeqMatrix](getSeqMatrix.md):         Calculate a matrix of pairwise Hamming distances 
for a sequence set.
+ [testSeqEqual](testSeqEqual.md):         Test sequences for equality.


Amino acid propertes
-------------------



+ [translateDNA](translateDNA.md):         Translate DNA sequences to amino acid sequences.
+ [aminoAcidProperties](aminoAcidProperties.md):  Calculate various physicochemical properties of amino acid 
sequences.
+ [countPatterns](countPatterns.md):        Count patterns in sequences.



General data manipulation
-------------------



+ [translateStrings](translateStrings.md):     Perform multiple string replacement operations.


References
-------------------


1. Vander Heiden JA, Yaari G, et al. pRESTO: a toolkit for processing 
high-throughput sequencing raw reads of lymphocyte receptor repertoires. 
Bioinformatics. 2014 30(13):1930-2.
1. Stern JNH, Yaari G, Vander Heiden JA, et al. B cells populating the multiple 
sclerosis brain mature in the draining cervical lymph nodes. 
Sci Transl Med. 2014 6(248):248ra107.
1. Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
peripheral blood IgE repertoires in patients with allergic rhinitis. 
J Allergy Clin Immunol. 2014 134(3):604-12.
1. Gupta NT, Vander Heiden JA, et al. Change-O: a toolkit for analyzing 
large-scale B cell immunoglobulin repertoire sequencing data.
Under review.





