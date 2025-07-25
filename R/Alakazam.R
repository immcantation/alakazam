#' @keywords internal
#' @aliases alakazam-package
"_PACKAGE"

# Alakazam package documentation and import directives

#' The Alakazam package
#' 
#' \code{alakazam} in a member of the Immcantation framework of tools and serves five main 
#' purposes:
#' \itemize{
#'   \item  Providing core functionality for other R packages in Immcantation. This
#'          includes common tasks such as file I/O, basic DNA sequence manipulation, and
#'          interacting with V(D)J segment and gene annotations.
#'   \item  Providing an R interface for interacting with the output of the pRESTO and 
#'          Change-O tool suites.
#'   \item  Performing clonal abundance and diversity analysis on lymphocyte repertoires.
#'   \item  Performing lineage reconstruction on clonal populations of immunoglobulin 
#'          (Ig) sequences.
#'   \item  Performing physicochemical property analyses of lymphocyte receptor sequences.
#' }
#' For additional details regarding the use of the \code{alakazam} package see the 
#' vignettes:\cr
#' \code{browseVignettes("alakazam")}
#' 
#' @section  File I/O:
#' \itemize{
#'   \item  \link{readChangeoDb}:        Input Change-O style files.
#'   \item  \link{writeChangeoDb}:       Output Change-O style files.
#' }
#' 
#' @section  Sequence cleaning:
#' \itemize{
#'   \item  \link{maskSeqEnds}:          Mask ragged ends.
#'   \item  \link{maskSeqGaps}:          Mask gap characters.
#'   \item  \link{collapseDuplicates}:   Remove duplicate sequences.
#' }
#' 
#' @section  Lineage reconstruction:
#' \itemize{
#'   \item  \link{makeChangeoClone}:     Clean sequences for lineage reconstruction.
#'   \item  \link{buildPhylipLineage}:   Perform lineage reconstruction of Ig sequences.
#' }
#' 
#' @section  Lineage topology analysis:
#' \itemize{
#'   \item  \link{tableEdges}:           Tabulate annotation relationships over edges.
#'   \item  \link{testEdges}:            Significance testing of annotation edges.
#'   \item  \link{testMRCA}:             Significance testing of MRCA annotations.
#'   \item  \link{summarizeSubtrees}:    Various summary statistics for subtrees.
#'   \item  \link{plotSubtrees}:         Plot distributions of summary statistics 
#'                                       for a population of trees.
#' }
#' 
#' @section  Diversity analysis:
#' \itemize{
#'   \item  \link{countClones}:          Calculate clonal abundance.
#'   \item  \link{estimateAbundance}:  	 Bootstrap clonal abundance curves.
#'   \item  \link{alphaDiversity}:  	 Generate clonal alpha diversity curves.
#'   \item  \link{plotAbundanceCurve}:   Plot clone size distribution as a rank-abundance 
#'   \item  \link{plotDiversityCurve}:   Plot clonal diversity curves.
#'   \item  \link{plotDiversityTest}:    Plot testing at given diversity hill indices. 
#' }
#' 
#' @section  Ig and TCR sequence annotation:
#' \itemize{
#'   \item  \link{countGenes}:           Calculate Ig and TCR allele, gene and family usage.
#'   \item  \link{extractVRegion}:       Extract CDRs and FWRs sub-sequences.
#'   \item  \link{getAllele}:            Get V(D)J allele names.
#'   \item  \link{getGene}:              Get V(D)J gene names.
#'   \item  \link{getFamily}:            Get V(D)J family names.
#'   \item  \link{junctionAlignment}: Junction alignment properties 
#' }
#' 
#' @section  Sequence distance calculation:
#' \itemize{
#'   \item  \link{seqDist}:        Calculate Hamming distance between two sequences.
#'   \item  \link{seqEqual}:       Test two sequences for equivalence.
#'   \item  \link{pairwiseDist}:   Calculate a matrix of pairwise Hamming distances for a 
#'                                 set of sequences.
#'   \item  \link{pairwiseEqual}:  Calculate a logical matrix of pairwise equivalence for a 
#'                                 set of sequences.
#' }
#' 
#' @section  Amino acid properties:
#' \itemize{
#'   \item  \link{translateDNA}:         Translate DNA sequences to amino acid sequences.
#'   \item  \link{aminoAcidProperties}:  Calculate various physicochemical properties of amino acid 
#'                                       sequences.
#'   \item  \link{countPatterns}:        Count patterns in sequences.
#'                                              
#' }
#' 
#' @name     alakazam
#' @references
#' \enumerate{
#'   \item  Vander Heiden JA, Yaari G, et al. pRESTO: a toolkit for processing 
#'            high-throughput sequencing raw reads of lymphocyte receptor repertoires. 
#'            Bioinformatics. 2014 30(13):1930-2.
#'   \item  Stern JNH, Yaari G, Vander Heiden JA, et al. B cells populating the multiple 
#'            sclerosis brain mature in the draining cervical lymph nodes. 
#'            Sci Transl Med. 2014 6(248):248ra107.
#'   \item  Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
#'            peripheral blood IgE repertoires in patients with allergic rhinitis. 
#'            J Allergy Clin Immunol. 2014 134(3):604-12.
#'   \item  Gupta NT, Vander Heiden JA, et al. Change-O: a toolkit for analyzing 
#'            large-scale B cell immunoglobulin repertoire sequencing data.
#'            Bioinformatics. 2015 Oct 15;31(20):3356-8.
#' }
#'
#' @import      ggplot2
#' @import      graphics
#' @import      methods
#' @import      utils
#' @importFrom  airr        read_rearrangement write_rearrangement
#' @importFrom	ape 		read.fastq read.tree di2multi reorder.phylo root ladderize
#' @importFrom  dplyr       do n desc %>%
#'                          bind_cols bind_rows combine arrange left_join
#'                          group_by ungroup
#'                          filter slice select 
#'                          mutate mutate_at 
#' 							one_of
#'							right_join rowwise
#'                          summarize summarize_at all_of
#'                          transmute rename
#' @importFrom  igraph      V E graph_from_data_frame as_data_frame as_edgelist 
#'                          make_graph make_directed_graph make_undirected_graph
#'                          vertex_attr set_vertex_attr
#'                          degree shortest_paths all_shortest_paths distances
#'                          graph_from_adjacency_matrix components groups
#' @importFrom  Matrix      sparseMatrix rowSums
#' @importFrom  progress    progress_bar
#' @importFrom  readr       read_delim read_tsv write_delim write_tsv cols
#' @importFrom  rlang       := sym syms enquo
#' @importFrom  scales      log2_trans log10_trans trans_breaks trans_format
#'                          math_format percent scientific pretty_breaks
#' @importFrom  seqinr      translate s2c
#' @importFrom  stats       na.omit setNames ecdf sd cor cov median mad
#'                          dbinom pbinom qbinom rbinom
#'                          dnorm pnorm qnorm rnorm
#'                          dmultinom rmultinom
#' @importFrom  stringi     stri_dup stri_flatten stri_join stri_length 
#'                          stri_count_boundaries stri_count_fixed 
#'                          stri_count_regex stri_extract_all_regex 
#'                          stri_extract_first_regex stri_replace_all_regex 
#'                          stri_replace_first_regex stri_split_fixed
#'                          stri_pad_left stri_pad_right
#'                          stri_detect_fixed stri_paste
#' @importFrom  tibble      tibble
#' @importFrom  tidyr       complete gather
#' @importFrom  Rcpp        evalCpp
#' @importFrom  Biostrings         BString extractAt
#' @importFrom  GenomicAlignments  explodeCigarOps explodeCigarOpLengths                       
#' @importFrom  IRanges            IRanges
#' @useDynLib   alakazam, .registration=TRUE
NULL

# Package loading actions
.onAttach <- function(libname, pkgname) {
    msg <- citation(pkgname)
    msg <-paste(c(format(msg,"citation")),collapse="\n\n")
    packageStartupMessage(msg)
}
