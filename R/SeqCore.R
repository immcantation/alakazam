# Common DNA, amino acid, and gene annotation operations for Alakazam
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2015.01.10


#### Constants ####

#' Default colors
#' 
#' Default color sets for nucleotide characters, Ig isotypes, and TCR chains.
#' 
#' @format  Named vectors with hexcode colors as values.
#' \itemize{
#'   \item \code{NUC_COLORS}  DNA character colors (A, C, G, T).
#'   \item \code{IG_COLORS}   Ig isotype colors (IgA, IgD, IgE, IgG, IgM, IgK, IgL).
#'   \item \code{TR_COLORS}   TCR chain colors (TRA, TRB, TRD, TRG).
#' }
#' 
#' @name   NUC_COLORS
#' @export
NUC_COLORS <- c("A"="#64F73F", 
                "C"="#FFB340", 
                "G"="#EB413C", 
                "T"="#3C88EE")

#' @rdname NUC_COLORS
#' @export
IG_COLORS <- c("IgA"="#377EB8", 
               "IgD"="#FF7F00", 
               "IgE"="#E41A1C", 
               "IgG"="#4DAF4A", 
               "IgM"="#984EA3",
               "IgK"="#E5C494",
               "IgL"="#FFD92F")

#' @rdname NUC_COLORS
#' @export
TR_COLORS <- c("TRA"="#CBD5E8", 
               "TRB"="#F4CAE4", 
               "TRD"="#FDCDAC", 
               "TRG"="#E6F5C9")


#' IUPAC ambiguous nucleotide translation
#'
#' A translation list for IUPAC ambiguous nucleotide characters and a corresponding vector
#' of matching DNA characters
#' 
#' @format  A list with single character codes as names and values containing character 
#'          vectors that containg the set of (A, C, G, T) characters corresponding to each 
#'          ambiguous character.
#' 
#' @name    IUPAC_NUC
#' @export
IUPAC_NUC <- list("A"="A", 
                  "C"="C", 
                  "G"="G", 
                  "T"="T",
                  "M"=c("A","C"), 
                  "R"=c("A","G"), 
                  "W"=c("A","T"), 
                  "S"=c("C","G"), 
                  "Y"=c("C","T"), 
                  "K"=c("G","T"), 
                  "V"=c("A","C","G"), 
                  "H"=c("A","C","T"), 
                  "D"=c("A","G","T"), 
                  "B"=c("C","G","T"),
                  "N"=c("A","C","G","T"))


#### Sequence distance functions ####

#' Build a nucleotide distance matrix
#'
#' \code{getNucMatrix} returns a Hamming distance matrix for IUPAC ambiguous
#' nucleotide characters with modifications for gap (-, .) and missing (?)
#' character values.
#' 
#' @param    gap  value to assign to gap (-, .) characters. 
#' @return   A \code{matrix} of nucleotide character distances with
#'           row and column names indicating the character pair. By default, 
#'           distances will be either O (equivalent), 1 (non-equivalent or missing), 
#'           or \code{Inf} (gap). 
#' 
#' @seealso  Creates nucleotide distance matrix for \code{\link{getSeqDistance}}.
#' @examples
#' # Set gap characters to Inf distance
#' # Distinguishes gaps from Ns
#' nuc_dist <- getNucMatrix()
#' 
#' # Set gap characters to 0 distance
#' # Makes gap characters equivalent to Ns
#' nuc_dist <- getNucMatrix(gap=0)
#' 
#' @export
getNucMatrix <- function(gap=Inf) {
    # Define Hamming distance matrix
    sub_mat <- diag(18)
    colnames(sub_mat) <- rownames(sub_mat) <- c(names(IUPAC_NUC), c("-", ".", "?"))
    for (i in 1:length(IUPAC_NUC)) {
        for (j in i:length(IUPAC_NUC)) {
            sub_mat[i, j] <- sub_mat[j, i] <- any(IUPAC_NUC[[i]] %in% IUPAC_NUC[[j]])
        }
    }
    
    # Add gap characters
    sub_mat[c(".", "-"), c(".", "-")] <- 1 
    sub_mat[c(".", "-"), 1:15] <- 1 - gap 
    sub_mat[1:15, c(".", "-")] <- 1 - gap
    
    return(1 - sub_mat)
}


#' Calculate distance between two sequences
#' 
#' \code{getSeqDistance} calculates the distance between two nucleotide
#' sequences.
#'
#' @param    seq1     character string containing a nucleotide sequence.
#' @param    seq2     character string containing a nucleotide sequence.
#' @param    nuc_mat  nucleotide character distance matrix. Defaults to
#'                    a Hamming distance matrix returned by \code{getNucMatrix}.
#'                    If gap characters (-, .) are assigned a value of \code{Inf}
#'                    in \code{nuc_mat} then continguous gaps of any run length,
#'                    which are not present in both sequences, will be counted as a 
#'                    distance of 1. Meaning, indels of any length will increase
#'                    the sequence distance by 1. Gap values other than \code{Inf}
#'                    will return a distance that does not consider indels as a 
#'                    special case.
#' @return   Distance between \code{seq1} and \code{seq2}. 
#' 
#' @seealso  Nucleotide distance matrix may be built with 
#'           \code{\link{getNucMatrix}}.
#' @examples
#' # Ungapped examples
#' getSeqDistance("ATGGC", "ATGGG")
#' getSeqDistance("ATGGC", "ATG??")
#' 
#' # Gaps will be treated as Ns with a gap=0 distance matrix
#' getSeqDistance("ATGGC", "AT--C", nuc_mat=getNucMatrix(gap=0))
#' 
#' # Gaps will be treated as universally non-matching characters with gap=1
#' getSeqDistance("ATGGC", "AT--C", nuc_mat=getNucMatrix(gap=1))
#' 
#' # Gaps of any length will be treated as single mismatches with a gap=Inf distance matrix
#' getSeqDistance("ATGGC", "AT--C", nuc_mat=getNucMatrix(gap=Inf))
#' 
#' # Gaps of equivalent run lengths are not counted as gaps
#' getSeqDistance("ATG-C", "ATG-C", nuc_mat=getNucMatrix(gap=Inf))
#'
#' # Overlapping runs of gap characters are counted as a single gap
#' getSeqDistance("ATG-C", "AT--C", nuc_mat=getNucMatrix(gap=Inf))
#' getSeqDistance("A-GGC", "AT--C", nuc_mat=getNucMatrix(gap=Inf))
#' getSeqDistance("AT--C", "AT--C", nuc_mat=getNucMatrix(gap=Inf))
#' 
#' # Discontiguous runs of gap characters each count as separate gaps
#' getSeqDistance("-TGGC", "AT--C", nuc_mat=getNucMatrix(gap=Inf))
#' 
#' @export
getSeqDistance <- function(seq1, seq2, nuc_mat=getNucMatrix(gap=Inf)) {
    # Convert string to character vector
    seq1 <- unlist(strsplit(seq1, ""))
    seq2 <- unlist(strsplit(seq2, ""))
    valid.idx <- (seq1 != "-" | seq2 != "-")
    seq1 <- seq1[valid.idx]
    seq2 <- seq2[valid.idx]
    # Calculate distance
    d <- sapply(1:length(seq1), function(x) { nuc_mat[seq1[x], seq2[x]] })
    indels <- sum(rle(d)$values == Inf)
    
    return(sum(d[is.finite(d)]) + indels)
}


#' Test nucleotide sequences for equality.
#' 
#' \code{testSeqEqual} checks if two nucleotide sequences are identical.
#'
#' @param    seq1    character string containing a nucleotide sequence.
#' @param    seq2    character string containing a nucleotide sequence.
#' @param    ignore  vector of characters to ignore when testing for equality.
#' @return   Returns \code{TRUE} if sequences are equal and \code{FALSE} if they are not.
#' 
#' @seealso  \code{\link{collapseDuplicates}}.
#' @examples
#' # Ignore gaps
#' testSeqEqual("ATG-C", "AT--C")
#' testSeqEqual("ATGGC", "ATGGN")
#' testSeqEqual("AT--T", "ATGGC")
#' 
#' # Ignore only Ns
#' testSeqEqual("ATG-C", "AT--C", ignore="N")
#' testSeqEqual("ATGGC", "ATGGN", ignore="N")
#' testSeqEqual("AT--T", "ATGGC", ignore="N")
#' 
#' @export
testSeqEqual <- function(seq1, seq2, ignore=c("N", "-", ".", "?")) {
    # Test that sequences lengths are equal
    if (nchar(seq1) != nchar(seq2)) {
        return(FALSE)
    }
    
    # Convert string to character vector
    x <- unlist(strsplit(seq1, ""))
    y <- unlist(strsplit(seq2, ""))
    
    # Determine non-ignored positions
    i <- !((x %in% ignore) | (y %in% ignore))
    
    return(all(x[i] == y[i]))
}


#### Sequence manipulation functions ####

#' Replace gap characters with Ns in nucleotide sequences
#' 
#' \code{maskSeqGaps} substitutes gap characters (-, .) with Ns in a
#' vector of nucleotide sequences.
#'
#' @param    seq         a character vector of nucleotide sequence strings.
#' @param    outer_only  if \code{TRUE} replace only continguous leading and trailing gaps;
#'                       if \code{FALSE} replace all gap characters.
#' @return   A modified \code{seq} vector with Ns in place of gap (-, .) characters.
#' 
#' @family   sequence manipulation functions
#' @seealso  Uses \code{\link{regex}} for replacement.
#'           
#' @examples
#' maskSeqGaps(c("ATG-C", "CC..C"))
#' maskSeqGaps("--ATG-C-")
#' maskSeqGaps("--ATG-C-", outer_only=TRUE)
#' 
#' @export
maskSeqGaps <- function(seq, outer_only=FALSE) {
    if (outer_only) {
        for (i in 1:length(seq)) {
            head_match <- attr(regexpr('^[-\\.]+', seq[i]), 'match.length')
            tail_match <- attr(regexpr('[-\\.]+$', seq[i]), 'match.length')
            if (head_match > 0) { 
                seq[i] <- gsub('^[-\\.]+', 
                                     paste(rep('N', head_match), collapse=''), 
                                     seq[i]) 
            }
            if (tail_match > 0) { 
                seq[i] <- gsub('[-\\.]+$', 
                                     paste(rep('N', tail_match), collapse=''), 
                                     seq[i]) 
            }
        }
    } else {
        seq <- gsub('[-\\.]', 'N', seq)
    }
    
    return(seq)
}


#' Replaces ragged leading and trailing edges of aligned sequences with Ns
#' 
#' \code{maskSeqEnds} takes a vector of nucleotide sequences, as character strings,
#' and replaces the leading and trailing characters with Ns to create a sequence
#' vector with uniformly masked outer sequence segments.
#' 
#' @param    seq       a character vector of nucleotide sequence strings.
#' @param    max_mask  the maximum number of characters to mask. If set to 0 then
#'                     no masking will be performed. If set to \code{NULL} then the upper 
#'                     masking bound will be automatically determined from the maximum 
#'                     number of observed leading or trailing Ns amongst all strings 
#'                     in \code{seq}. 
#' @param    trim      if \code{TRUE} leading and trailing characters will be cut rather 
#'                     than masked with Ns.
#' @return   A modified \code{seq} vector with masked (or optionally trimmed) sequences.
#' 
#' @family   sequence manipulation functions
#' @seealso  Uses \code{\link{regex}} and \code{\link{substr}} for replacement.
#' @examples
#' # Default behavior uniformly masks ragged ends
#' seq <- c("CCCCTGGG", "NAACTGGN", "NNNCTGNN")
#' maskSeqEnds(seq)
#'
#' # Does nothing
#' maskSeqEnds(seq, max_mask=0)
#' 
#' # Cut ragged sequence ends
#' maskSeqEnds(seq, trim=TRUE)
#'
#' # Set max_mask to limit extent of masking and trimming
#' maskSeqEnds(seq, max_mask=1)
#' maskSeqEnds(seq, max_mask=1, trim=TRUE)
#' 
#' @export
maskSeqEnds <- function(seq, max_mask=NULL, trim=FALSE) {
    # Find length of leading and trailing Ns
    left_lengths <- attr(regexpr('(^N*)', seq, perl=T), 'capture.length')
    right_lengths <- attr(regexpr('(N*$)', seq, perl=T), 'capture.length')
    
    # Mask to minimal inner sequence length
    left_mask <- min(max(left_lengths[, 1]), max_mask)
    right_mask <- min(max(right_lengths[, 1]), max_mask)
    seq_lengths <- nchar(seq)
    if (trim) {
        seq <- substr(seq, left_mask + 1, seq_lengths - right_mask)
    } else {
        substr(seq, 0, left_mask) <- paste(rep('N', left_mask), collapse='')
        substr(seq, seq_lengths - right_mask + 1, seq_lengths + 1) <- paste(rep('N', right_mask), collapse='')
    }
    
    return(seq)
}


#' Remove duplicate nucleotide sequences and combine annotations
#' 
#' \code{collapseDuplicates} identifies duplicate nucleotides sequences in a 
#' data.frame, allowing for ambiguous nucleotide characters, removes the
#' duplicate entries and combines any associated annotations (data.frame columns).
#'
#' @param    data         data.frame containing Change-O columns. The data.frame frame
#'                        must contain, at a minimum, a unique identifier column 
#'                        and a column containg a character vector of nucleotide sequences.
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing nucleotide sequences.
#' @param    text_fields  character vector of textual columns to collapse. The textual 
#'                        annotations of duplicate sequences will be merged into a single 
#'                        string with each unique value alphabetized and delimited by a 
#'                        "/" character.
#' @param    num_fields   a vector of numeric columns to collapse. The numeric annotations
#'                        of duplicate sequences will be summed. 
#' @param    seq_fields   a vector of nucletoide sequence columns to collapse. The sequence 
#'                        with the fewest Ns will be retained. This is distinct from the 
#'                        \code{seq} parameter which is used to determine duplicates.
#' @param    ignore       vector of characters to ignore when testing for equality.
#' @param    verbose      if \code{TRUE} report the number input, discarded and output 
#'                        sequences; if \code{FALSE} process sequences silently.
#' @return   A modified data.frame with duplicate sequences removed and annotation fields 
#'           collapsed. Columns that are not specified in either \code{text_fields} or
#'           \code{num_fields} will be retained, but the value will be chosen from a random
#'           sequences amongst all sequences in a cluster of duplicates.  Sequences that
#'           could not be unambiguously assigned to a single duplicate cluster are discarded,
#'           along with their annotations.
#' 
#' @family   sequence manipulation functions
#' @seealso  \code{\link{testSeqEqual}}.
#' @examples
#' # Example Change-O data.frame
#' df <- data.frame(SEQUENCE_ID=LETTERS[1:4],
#'                  SEQUENCE_GAP=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
#'                  TYPE=c("IgM", "IgG", "IgG", "IgA"),
#'                  SAMPLE=c("S1", "S1", "S2", "S2"),
#'                  COUNT=1:4,
#'                  stringsAsFactors=FALSE)
#' 
#' # Annotations are not parsed if neither text_fields nor num_fields is specified
#' # The retained sequence annotations will be random
#' collapseDuplicates(df, verbose=T)
#'
#' # Unique text_fields annotations are combined into a single string with "/"
#' # num_fields annotations are summed
#' # Unambiguous duplicates are discarded
#' collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", verbose=T)
#' 
#' # Masking ragged ends may impact duplicate removal
#' df$SEQUENCE_GAP <- maskSeqEnds(df$SEQUENCE_GAP)
#' collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", verbose=T)
#'
#' @export
collapseDuplicates <- function(data, id="SEQUENCE_ID", seq="SEQUENCE_GAP",
                               text_fields=NULL, num_fields=NULL, seq_fields=NULL,
                               ignore=c("N", "-", ".", "?"), verbose=FALSE) {
    # Verify column classes and exit if they are incorrect
    if (!all(sapply(subset(data, select=text_fields), is.character))) {
        stop("All text_fields columns must be of type 'character'")
    }
    if (!all(sapply(subset(data, select=num_fields), is.numeric))) {
        stop("All num_fields columns must be of type 'numeric'")
    }
    if (!all(sapply(subset(data, select=seq_fields), is.character))) {
        stop("All seq_fields columns must be of type 'character'")
    }
    seq_len <- nchar(data[, seq])
    if (any(seq_len != seq_len[1])) {
        warning("All sequences are not the same length")
    }
    
    # Define verbose reporting function
    printVerbose <- function(n_total, n_unique, n_discard) {
        cat(" FUNCTION> collapseDuplicates\n", sep="")
        cat("    TOTAL> ", n_total, "\n", sep="")
        cat("   UNIQUE> ", n_unique, "\n", sep="")
        cat("COLLAPSED> ", n_total - n_unique - n_discard, "\n", sep="")
        cat("DISCARDED> ", n_discard, "\n", sep="")
        cat("\n")
    }
    
    # Return input if there are no sequences to collapse
    nseq <- nrow(data)
    if (nseq <= 1) { 
        if (verbose) { printVerbose(nseq, 1, 0) }
        return(data)
    }
    
    # Build distance matrix
    d_mat <- matrix(TRUE, nseq, nseq, dimnames=list(data[, id], data[, id]))
    for (i in 1:(nseq - 1)) {
        for (j in (i + 1):nseq) {
            d_mat[i, j] <- d_mat[j, i] <- testSeqEqual(data[i, seq], data[j, seq], ignore)
        }
    }
    
    # Return input if no sequences are equal
    if (!any(d_mat[lower.tri(d_mat, diag=F)])) {
        if (verbose) { printVerbose(nseq, nseq, 0) }
        return(data)
    }        
    
    # Find sequences that will cluster ambiguously
    ambig_rows <- numeric()
    for (i in 1:nseq) {
        idx <- which(d_mat[i, ])
        tmp_mat <- d_mat[idx, idx]
        if (!all(tmp_mat)) { 
            ambig_rows <- append(ambig_rows, i) 
        }
    }
    discard_count <- length(ambig_rows)

    # Return single sequence if all sequence belong to ambiguous clusters
    if (discard_count == nrow(d_mat)) {
        unambig_len <- nchar(gsub("N", "", data[, seq]))
        if (verbose) { printVerbose(nseq, 0, discard_count - 1) }
        return(data[which.max(unambig_len), ])
    }
    
    # Exclude ambiguous sequences from clustering
    if (discard_count > 0) {
        d_mat <- d_mat[-ambig_rows, -ambig_rows]
    }
        
    # Cluster remaining sequences into unique and duplicate sets
    dup_taxa <-  list()
    uniq_taxa <- character()
    done_taxa <- character()
    taxa_names <- rownames(d_mat)
    for (taxa in taxa_names) {
        # Skip taxa if previously assigned to a cluster
        if (taxa %in% done_taxa) { next }
            
        # Find all zero distance taxa
        idx <- which(d_mat[taxa, ])
        if (length(idx) == 1) {
            # Assign unique sequences to unique vector
            uniq_taxa <- append(uniq_taxa, taxa_names[idx])
        } else if (length(idx) > 1) {
            # Assign clusters of duplicates to duplicate list            
            dup_taxa <- c(dup_taxa, list(taxa_names[idx]))    
        } else {
            # Report error (should never occur)
            stop("Error in distance matrix of collapseDuplicates")
        }
        # Update vector of clustered taxa
        done_taxa <- c(done_taxa, taxa_names[idx])
    }
    
    # Collapse duplicate sets and append entries to unique data.frame
    unique_list <- list(data[data[, id] %in% uniq_taxa, ])
    for (taxa in dup_taxa) {
        # Define row indices of identical sequences
        idx <- which(data[, id] %in% taxa)
        tmp_df <- data[idx[1], ]
        
        if (length(idx) > 1) {
            # Define set of text fields for row
            for (f in text_fields) {
                f_set <- na.omit(data[idx, f])
                if (length(f_set) > 0) {
                    f_set <- unlist(strsplit(f_set, '/'))
                    f_set <- sort(unique(f_set))
                    f_val <- paste(f_set, collapse='/')
                } else {
                    f_val <- NA
                }
                tmp_df[, f] <- f_val
            }
            
            # Sum numeric fields
            for (f in num_fields) {
                f_set <- na.omit(data[idx, f])
                if (length(f_set) > 0) { 
                    f_val <- sum(f_set) 
                } else { 
                    f_val <- NA 
                }
                tmp_df[, f] <- f_val
            }
            
            # Select sequence fields with fewest Ns
            for (f in seq_fields) {
                f_set <- na.omit(data[idx, f])
                if (length(f_set) > 0) {
                    f_len <- nchar(gsub("N", "", f_set))
                    f_val <- f_set[which.max(f_len)]
                } else {
                    f_val <- NA
                }
                tmp_df[, f] <- f_val
            }
            
            # Assign id and sequence with least number of Ns
            seq_set <- data[idx, c(id, seq)]
            unambig_len <- nchar(gsub("N", "", seq_set[, seq]))
            tmp_df[, c(id, seq)] <- seq_set[which.max(unambig_len), c(id, seq)]
        }
        
        # Add row to unique list
        unique_list <- c(unique_list, list(tmp_df))
    }
    
    # Combine all rows into unique data.frame
    unique_df <- plyr::rbind.fill(unique_list)
    
    if (verbose) { printVerbose(nseq, nrow(unique_df), discard_count) }
    return(unique_df)
}


#' Extracts FWRs and CDRs from IMGT-gapped sequences
#' 
#' \code{getVRegion} extracts the framework and complementarity determining regions of the V-segment 
#' for IMGT-gapped immunoglobulin (Ig) nucleotide sequences from  according to the IMGT numbering scheme.
#'
#' @param     sequences  character vector of IMGT gapped sequences.
#' @param     region     string defining the region of the V segment to extract. Must be one of 
#'                       \code{c("FWR1", "CDR1", "FWR2", "CDR2" ,"FWR3")}
#' @return    A character vector of the extracted sub-sequences.
#' 
#' @seealso   \code{\link{substr}}.
#' @references
#'   \url{http://imgt.org}
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' clone <- subset(df, CLONE == 164)
#'
#' # Get regions
#' getVRegion(clone$SEQUENCE_GAP, "FWR1")
#' getVRegion(clone$SEQUENCE_GAP, "CDR1")
#' getVRegion(clone$SEQUENCE_GAP, "FWR2")
#' getVRegion(clone$SEQUENCE_GAP, "CDR2")
#' getVRegion(clone$SEQUENCE_GAP, "FWR3")
#'
#' @export
getVRegion <- function(sequences, region) {
    # Define IMGT region boundaries
    imgt <- list("FWR1"=c(1, 78),
                 "CDR1"=c(79, 114),
                 "FWR2"=c(115, 165),
                 "CDR2"=c(166, 195),
                 "FWR3"=c(196, 312))
    # Check region argument    
    if (!(region %in% names(imgt))) {
        stop("Incorrect region value. Must be one of c('FWR1', 'CDR1', 'FWR2', 'CDR2' ,'FWR3').")
    }
    
    return(substr(sequences, imgt[[region]][1], imgt[[region]][2]))
}

  
#### Gene annotation functions ####

#' Get Ig segment allele, gene and family names
#' 
#' \code{getSegment} performs generic matching of delimited segment calls with a custom regular 
#' expression. \code{getAllele}, \code{getGene} and \code{getFamily} extract the allele, gene and 
#' family names, respectively, from a character vector of immunoglobulin (Ig) segment allele calls 
#' in IMGT format.
#'
#' @param     segment_call    character vector containing segment calls delimited by commas.
#' @param     segment_regex   string defining the segment match regular expression.
#' @param     first           if \code{TRUE} return only the first call in \code{segment_call};
#'                            if \code{FALSE} return all calls delimited by commas.
#' @param     collapse        if \code{TRUE} check for duplicates and return only unique segment
#'                            assignments; if \code{FALSE} return all assignments (faster). 
#'                            Has no effect if \code{first=TRUE}.
#' @param     sep             character defining both the input and output segment call delimiter.
#' @return    A character vector containing allele, gene or family names
#' 
#' @seealso   Uses \code{\link{str_extract}}.
#' @references
#'   \url{http://imgt.org}
#' @examples
#' kappa_call <- c("Homsap IGKV1-39*01 F,Homsap IGKV1D-39*01 F", "Homsap IGKJ5*01 F")
#'
#' getAllele(kappa_call)
#' getAllele(kappa_call, first=FALSE)
#' 
#' getGene(kappa_call)
#' getGene(kappa_call, first=FALSE)
#' 
#' getFamily(kappa_call)
#' getFamily(kappa_call, first=FALSE)
#' getFamily(kappa_call, first=FALSE, collapse=TRUE)
#'
#' @export
getSegment <- function(segment_call, segment_regex, first=TRUE, collapse=TRUE, sep=",") {
    # Define boundaries of individual segment calls
    edge_regex <- if (first) { ".*" } else { paste0("[^", sep, "]*") }
    
    # Extract calls
    r <- gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), "\\1", 
              segment_call, perl=T)
    
    # Collapse to unique set if required
    if (!first & collapse) {
        r <- sapply(strsplit(r, sep), function(x) paste(unique(x), collapse=sep))
    }
    
    return(r)
}

#' @rdname getSegment
#' @export
getAllele <- function(segment_call, first=TRUE, collapse=TRUE, sep=",") {
    allele_regex <- '((IG[HLK]|TR[ABGD])[VDJ]\\d+[-/\\w]*[-\\*][\\.\\w]+)'
    r <- getSegment(segment_call, allele_regex, first=first, collapse=collapse, sep=sep)
    
    return(r)
}

#' @rdname getSegment
#' @export
getGene <- function(segment_call, first=TRUE, collapse=TRUE, sep=",") {
    gene_regex <- '((IG[HLK]|TR[ABGD])[VDJ]\\d+[-/\\w]*)'
    r <- getSegment(segment_call, gene_regex, first=first, collapse=collapse, sep=sep)
    
    return(r)
}


#' @rdname getSegment
#' @export
getFamily <- function(segment_call, first=TRUE, collapse=TRUE, sep=",") {
    family_regex <- '((IG[HLK]|TR[ABGD])[VDJ]\\d+)'
    r <- getSegment(segment_call, family_regex, first=first, collapse=collapse, sep=sep)
    
    return(r)
}




