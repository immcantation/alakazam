Alakazam
-------------------------------------------------------------------------------

Alakazam is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework for Adaptive Immune Receptor Repertoire sequencing 
(AIRR-seq) and provides a set of tools to investigate lymphocyte 
receptor clonal lineages, diversity, gene usage, and other repertoire level 
properties, with a focus on high-throughput immunoglobulin (Ig) sequencing.

Alakazam serves five main purposes:

1. Providing core functionality for other R packages in the Immcantation 
   framework. This includes common tasks such as file I/O, basic DNA sequence 
   manipulation, and interacting with V(D)J segment and gene annotations.
2. Providing an R interface for interacting with the output of the 
   [pRESTO](http://presto.readthedocs.io) and 
   [Change-O](http://changeo.readthedocs.io) tool suites.
3. Performing lineage reconstruction on clonal populations of Ig sequences 
   and analyzing the topology of the resultant lineage trees. 
4. Performing clonal abundance and diversity analysis on lymphocyte 
   repertoires.
5. Performing physicochemical property analyses of lymphocyte receptor 
   sequences.


Contact
-------------------------------------------------------------------------------

For help and questions please contact the [Immcantation Group](mailto:immcantation@googlegroups.com)
or use the [issue tracker](https://bitbucket.org/kleinstein/alakazam/issues?status=new&status=open).


# Dependencies

**Depends:** ggplot2  
**Imports:** ape, dplyr, graphics, grid, igraph, lazyeval, Matrix, methods, progress, Rcpp, readr, rlang, scales, seqinr, stats, stringi, tibble, tidyr, utils  
**Suggests:** knitr, rmarkdown, testthat  
**Extends:** FALSE


# Authors

[Immcantation](mailto:immcantation@gmail.com) (cre)  
[Jason Vander Heiden](mailto:jason.vanderheiden@yale.edu) (aut)  
[Namita Gupta](mailto:namita.gupta@yale.edu) (aut)  
[Susanna Marquez](mailto:susanna.marquez@yale.edu) (ctb)  
[Daniel Gadala-Maria](mailto:daniel.gadala-maria@yale.edu) (ctb)  
[Ruoyi Jiang](mailto:ruoyi.jiang@yale.edu) (ctb)  
[Nima Nouri](mailto:nima.nouri@yale.edu) (ctb)  
[Kenneth Hoehn](mailto:kenneth.hoehn@yale.edu) (ctb)  
[Julian Zhou](mailto:julian.zhou@bulldogs.yale.edu) (ctb)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


# Citing


To cite the alakazam package in publications, please use:

Gupta N, Vander Heiden J, Uduman M, Gadala-Maria D, Yaari G, Kleinstein S (2015). “Change-O: a toolkit for
analyzing large-scale B cell immunoglobulin repertoire sequencing data.” _Bioinformatics_, 1-3. doi:
10.1093/bioinformatics/btv359 (URL: https://doi.org/10.1093/bioinformatics/btv359).

To cite the Ig-specific lineage reconstruction and diversity methods, please use:

Stern J, Yaari G, Vander Heiden J, Church G, Donahue W, Hintzen R, Huttner A, Laman J, Nagra R, Nylander A, Pitt
D, Ramanan S, Siddiqui B, Vigneault F, Kleinstein S, Hafler D, O'Connor K (2014). “B cells populating the
multiple sclerosis brain mature in the draining cervical lymph nodes.” _Science Translational Medicine_,
*6*(248), 248ra107. doi: 10.1126/scitranslmed.3008879 (URL: https://doi.org/10.1126/scitranslmed.3008879).


