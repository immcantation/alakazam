# Project documentation for Alakazam
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2015.04.01


#' The alakazam package
#' 
#' \code{alakazam} in a member of the Change-O suite of tools and serves three main 
#' purposes:
#' \itemize{
#'   \item  Providing core functionality for other R packages in the Change-O suite. This
#'          includes common tasks such as file I/O and basic DNA sequence manipulation.
#'   \item  Providing an R interface for interacting with the output of the pRESTO tool suite.
#'   \item  Providing tools to perform lineage reconstruction and diversity analysis on clonal 
#'          populations of immunoglobulin (Ig) sequences. 
#' }
#' For additional details regarding the use of the \code{alakazam} package see the 
#' vignettes:\cr
#' \code{browseVignettes("alakazam")}
#' 
#' @section  File I/O:
#' \itemize{
#'   \item  \code{\link{readChangeoDb}}:  Input Change-O style files.
#'   \item  \code{\link{writeChangeoDb}}:  Output Change-O style files.
#' }
#' 
#' @section  Sequence cleaning:
#' \itemize{
#'   \item  \code{\link{maskSeqEnds}}:         Mask ragged ends.
#'   \item  \code{\link{maskSeqGaps}}:         Mask gap characters.
#'   \item  \code{\link{collapseDuplicates}}:  Remove duplicate sequences.
#' }
#' 
#' @section  Lineage reconstruction:
#' \itemize{
#'   \item  \code{\link{makeChangeoClone}}:    Clean sequences for lineage reconstruction.
#'   \item  \code{\link{buildPhylipLineage}}:  Perform lineage reconstruction of Ig sequences.
#' }
#' 
#' @section  Diversity analysis:
#' \itemize{
#'   \item  \code{\link{bootstrapDiversity}}:  Generate clonal diversity curves.
#'   \item  \code{\link{plotDiversityCurve}}:  Plot clonal diversity curves.
#'   \item  \code{\link{testDiversity}}:       Test significance of clonal diversity scores.
#' }
#' 
#' @section  Ig sequence annotation:
#' \itemize{
#'   \item  \code{\link{extractVRegion}}:      Extract CDRs and FWRs sub-sequences.
#'   \item  \code{\link{getAllele}}:           Get V(D)J allele names.
#'   \item  \code{\link{getGene}}:             Get V(D)J gene names.
#'   \item  \code{\link{getFamily}}:           Get V(D)J family names.
#' }
#' 
#' @section  Sequence distance calculation:
#' \itemize{
#'   \item  \code{\link{getSeqDistance}}:      Calculate Hamming distance between sequences.
#'   \item  \code{\link{testSeqEqual}}:        Test sequences for equality.
#' }
#' 
#' @section  General data manipulation:
#' \itemize{
#'   \item  \code{\link{translateStrings}}:    Perform multiple string replacement operations.
#' } 
#' 
#' @seealso
#' The Change-O suite of tools includes three separate R packages: \code{\link[alakazam]{alakazam}}, 
#' \code{\link[tigger]{tigger}}, and \code{\link[shm]{shm}}.
#' 
#' @name     alakazam
#' @docType  package
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
#'            Under review.
#' }
#' 
#' @import   ggplot2
#' @import   igraph
#' @import   plyr
#' @import   scales
NULL
