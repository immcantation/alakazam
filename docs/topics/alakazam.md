# The Alakazam package

Description
--------------------

`alakazam` in a member of the Immcantation framework of tools and serves five main 
purposes:

+ Providing core functionality for other R packages in Immcantation. This
includes common tasks such as file I/O, basic DNA sequence manipulation, and
interacting with V(D)J segment and gene annotations.
+ Providing an R interface for interacting with the output of the pRESTO and 
Change-O tool suites.
+ Performing clonal abundance and diversity analysis on lymphocyte repertoires.
+ Performing lineage reconstruction on clonal populations of immunoglobulin 
(Ig) sequences.
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



Lineage topology analysis
-------------------



+ [tableEdges](tableEdges.md):           Tabulate annotation relationships over edges.
+ [testEdges](testEdges.md):            Significance testing of annotation edges.
+ [testMRCA](testMRCA.md):             Significance testing of MRCA annotations.
+ [summarizeSubtrees](summarizeSubtrees.md):    Various summary statistics for subtrees.
+ [plotSubtrees](plotSubtrees.md):         Plot distributions of summary statistics 
for a population of trees.



Diversity analysis
-------------------



+ [countClones](countClones.md):          Calculate clonal abundance.
+ [estimateAbundance](estimateAbundance.md):  	 Bootstrap clonal abundance curves.
+ [alphaDiversity](alphaDiversity.md):  	 Generate clonal alpha diversity curves.
+ [plotAbundanceCurve](plotAbundanceCurve.md):   Plot clone size distribution as a rank-abundance 
+ [plotDiversityCurve](plotDiversityCurve.md):   Plot clonal diversity curves.
+ [plotDiversityTest](plotDiversityTest.md):    Plot testing at given diversity hill indices. 



Ig and TCR sequence annotation
-------------------



+ [countGenes](countGenes.md):           Calculate Ig and TCR allele, gene and family usage.
+ [extractVRegion](extractVRegion.md):       Extract CDRs and FWRs sub-sequences.
+ [getAllele](getSegment.md):            Get V(D)J allele names.
+ [getGene](getSegment.md):              Get V(D)J gene names.
+ [getFamily](getSegment.md):            Get V(D)J family names.
+ [junctionAlignment](junctionAlignment.md): Junction alignment properties 



Sequence distance calculation
-------------------



+ [seqDist](seqDist.md):        Calculate Hamming distance between two sequences.
+ [seqEqual](seqEqual.md):       Test two sequences for equivalence.
+ [pairwiseDist](pairwiseDist.md):   Calculate a matrix of pairwise Hamming distances for a 
set of sequences.
+ [pairwiseEqual](pairwiseEqual.md):  Calculate a logical matrix of pairwise equivalence for a 
set of sequences.



Amino acid properties
-------------------



+ [translateDNA](translateDNA.md):         Translate DNA sequences to amino acid sequences.
+ [aminoAcidProperties](aminoAcidProperties.md):  Calculate various physicochemical properties of amino acid 
sequences.
+ [countPatterns](countPatterns.md):        Count patterns in sequences.




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
Bioinformatics. 2015 Oct 15;31(20):3356-8.





See also
-------------------

Useful links:

+  [https://alakazam.readthedocs.io/](https://alakazam.readthedocs.io/)
+  Report bugs at [https://github.com/immcantation/alakazam/issues](https://github.com/immcantation/alakazam/issues)





Author
-------------------

**Maintainer**: Susanna Marquez [susanna.marquez@yale.edu](susanna.marquez@yale.edu)

Authors:

+  Namita Gupta [namita.gupta@yale.edu](namita.gupta@yale.edu)
+  Nima Nouri [nima.nouri@yale.edu](nima.nouri@yale.edu)
+  Ruoyi Jiang [ruoyi.jiang@yale.edu](ruoyi.jiang@yale.edu)
+  Julian Zhou [julian.zhou@bulldogs.yale.edu](julian.zhou@bulldogs.yale.edu)
+  Kenneth Hoehn [kenneth.hoehn@yale.edu](kenneth.hoehn@yale.edu)
+  Jason Vander Heiden [jason.vanderheiden@gmail.com](jason.vanderheiden@gmail.com)
+  Steven Kleinstein [steven.kleinstein@yale.edu](steven.kleinstein@yale.edu) [copyright holder]


Other contributors:

+  Daniel Gadala-Maria [daniel.gadala-maria@yale.edu](daniel.gadala-maria@yale.edu) [contributor]
+  Edel Aron [edel.aron@yale.edu](edel.aron@yale.edu) [contributor]
+  Cole Jensen [cole.jensen@yale.edu](cole.jensen@yale.edu) [contributor]
+  Gisela Gabernet [gisela.gabernet@yale.edu](gisela.gabernet@yale.edu) [contributor]





