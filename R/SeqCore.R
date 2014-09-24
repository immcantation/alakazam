# Common DNA, amino acid, and gene annotation operations for Alakazam
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2014.9.24


#### Constants ####

#' Default colors
NUC_COLORS <- c("A"="#64F73F", 
                "C"="#FFB340", 
                "G"="#EB413C", 
                "T"="#3C88EE",
                "U"="#3C88EE")
IG_COLORS <- c("IgA"="#377EB8", 
               "IgD"="#FF7F00", 
               "IgE"="#E41A1C", 
               "IgG"="#4DAF4A", 
               "IgM"="#984EA3",
               "IgK"="#E5C494",
               "IgL"="#FFD92F")
TR_COLORS <- c("TRA"="#CBD5E8", 
               "TRB"="#F4CAE4", 
               "TRD"="#FDCDAC", 
               "TRG"="#E6F5C9")

#' IUPAC nucleotide translation list
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


#### Sequence manipulation functions ####

#' Replace gap characters with Ns in nucleotide sequences
#' 
#' \code{maskSeqGaps} substitutes gap characters (-, .) with Ns in a
#' vector of nucleotide sequences.
#'
#' @param    seq         a character vector of nucleotide sequence strings.
#' @param    outer_only  if TRUE replace only continguous leading and trailing gaps;
#'                       if FALSE replace all gap characters.
#' @return   A modified \code{seq} vector with Ns in place of gap (-, .) characters.
#' 
#' @family   sequence manipulation functions
#' @seealso  Called by \code{\link{prepLineageClone}}. 
#'           Uses \code{\link{regex}} for replacement.
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
#' @seealso  Called by \code{\link{prepLineageClone}}. 
#'           Uses \code{\link{regex}} and \code{\link{substr}} for replacement.
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
#' @param    nuc_mat      nucleotide character distance matrix.
#' @return   A modified data.frame with duplicate sequences removed and annotation fields 
#'           collapsed. Columns that are not specified in either \code{text_fields} or
#'           \code{num_fields} will be retained, but the value will be chosen from a random
#'           sequences amongst all sequences in a cluster of duplicates.  Sequences that
#'           could not be unambiguously assigned to a single duplicate cluster are discarded,
#'           along with their annotations.
#' 
#' @family   sequence manipulation functions
#' @seealso  Called by \code{\link{prepLineageClone}}. 
#'           Distance is determined using \code{\link{getSeqDistance}} and 
#'           \code{\link{getNucMatrix}}.
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
#' collapseDuplicates(df)
#'
#' # Unique text_fields annotations are combined into a single string with "/"
#' # num_fields annotations are summed
#' # Unambiguous duplicates are discarded
#' collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT")
#' 
#' # Masking ragged ends may impact duplicate removal
#' df$SEQUENCE_GAP <- maskSeqEnds(df$SEQUENCE_GAP)
#' collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT")
#'
#' @export
collapseDuplicates <- function(data, id="SEQUENCE_ID", seq="SEQUENCE_GAP",
                               text_fields=NULL, num_fields=NULL, 
                               nuc_mat=getNucMatrix(gap=0)) {
    
    # >>> REMOVE rbind calls for speed
    # >>> Remove or cleanup progress messages
    cat('FUNCTION> collapseDuplicates\n')
    
    # Return input if there are no sequences to collapse
    nseq <- nrow(data)
    if (nseq <= 1) { 
        cat('INFO>', nseq, 'sequences total\n')
        cat('INFO> 0 sequences collapsed\n')
        cat('INFO> 0 sequences discarded\n')
        cat('\n')
        return(data)
    }
    
    # Build distance matrix
    d_mat <- matrix(0, nseq, nseq, 
                    dimnames=list(data[, id], data[, id]))
    for (i in 1:(nseq - 1)) {
        for (j in (i + 1):nseq) {
            d_mat[i, j] <- d_mat[j, i] <- getSeqDistance(data[i, seq], 
                                                         data[j, seq], 
                                                         nuc_mat)
        }
    }
    
    # Return input if no sequences have zero distance
    if (all(d_mat[lower.tri(d_mat, diag=F)] != 0)) {
        cat('INFO>', nseq, 'sequences total\n')
        cat('INFO> 0 sequences collapsed\n')
        cat('INFO> 0 sequences discarded\n')
        cat('\n')
        return(data)
    }        
    
    # Find sequences that will cluster ambiguously
    ambig_rows <- numeric()
    for (i in 1:nseq) {
        idx <- which(d_mat[i, ] == 0)
        tmp_mat <- d_mat[idx, idx]
        if (any(tmp_mat != 0)) { 
            ambig_rows <- append(ambig_rows, i) 
        }
    }
    
    # Exclude ambiguous sequences from clustering
    discard_count <- length(ambig_rows)
    discard_ids <- rownames(d_mat)[ambig_rows]
    if (discard_count > 0) {
        d_mat <- d_mat[-ambig_rows, -ambig_rows]
    }
    
    # Cluster remaining sequences into unique and duplicate sets
    dup_taxa <-  list()
    uniq_taxa <- character()
    done_taxa <- character()
    taxa_names <- rownames(d_mat)
    for (i in 1:nrow(d_mat)) {
        # Skip taxa if previously assigned to a cluster
        if (taxa_names[i] %in% done_taxa) { 
            next 
        }
        
        # Find all zero distance taxa
        idx <- which(d_mat[i, ] == 0)
        if (length(idx) == 1) {
            # Assign unique sequences to unique vector         
            uniq_taxa <- append(uniq_taxa, taxa_names[idx])
        } else if (length(idx) > 1) {
            # Assign clusters of duplicates to duplicate list            
            dup_taxa <- append(dup_taxa, list(taxa_names[idx]))
        } else {
            # Report error (should never occur)
            stop('Error in distance matrix of collapseDuplicates')
        }
        # Update vector of clustered taxa
        done_taxa <- append(done_taxa, taxa_names[idx])
    }
    
    # Get data.frame of unique sequences
    unique_df <- data[data[, id] %in% uniq_taxa, ]
    
    # Collapse duplicate sets and append entries to unique data.frame
    for (taxa in dup_taxa) {
        # Define row indices of identical sequences
        idx <- which(data[, id] %in% taxa)
        tmp_df <- data[idx[1], ]
        
        if (length(idx) > 1) {
            # Define set of text fields for row
            for (f in text_fields) {
                f_set <- na.omit(data[idx, f])
                f_set <- unlist(strsplit(f_set, '/'))
                f_set <- sort(unique(f_set))
                tmp_df[f] <- paste(f_set, collapse='/')
            }
            
            # Sum numeric fields
            for (f in num_fields) {
                f_set <- data[idx, f]
                tmp_df[f] <- sum(f_set, na.rm=T)
            }
            
            #Assign sequence with least number of N characters
            seq_set <- unique(data[idx, seq])
            unambig_len <- nchar(gsub("N", "", seq_set))
            tmp_df[, seq] <- seq_set[which.max(unambig_len)]
        }
        
        # Add row to unique data.frame
        unique_df <- plyr::rbind.fill(unique_df, tmp_df)
    }
    
    cat('INFO>', nseq, 'sequences total\n')
    cat('INFO>', nseq - nrow(unique_df) - discard_count, 'sequences collapsed\n')
    if (discard_count > 0) {
        cat('INFO> ', discard_count, ' sequences discarded (', paste(discard_ids, collapse=','), ')\n', sep='')
    } else {
        cat('INFO> 0 sequences discarded\n')
    }
    cat('\n')
    
    return(unique_df)
}


#### Gene annotation functions ####

#' Get Ig segment allele, gene and family names
#' 
#' \code{getAllele}, \code{getGene} and \code{getFamily} extract the allele, gene and family names,
#' respectively, from a character vector of immunoglobulin (Ig) segment allele calls in IMGT format.
#'
#' @param     segment_call    character vector containing segment calls delimited by commas.
#' @param     first           if TRUE return only the first call in \code{segment_call};
#'                            if FALSE return all calls delimited by commas.
#' @return    a character vector containing allele, gene or family names
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
#'
#' @export
getAllele <- function(segment_call, first=TRUE) {
    allele_regex <- '(IG[HLK][VDJ]\\d+[-/\\w]*[-\\*][\\.\\w]+)'
    
    if (first) {
        r <- stringr::str_extract(segment_call, stringr::perl(allele_regex))
    } else {
        r <- stringr::str_extract_all(segment_call, stringr::perl(allele_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}

#' @rdname getAllele
#' @export
getGene <- function(segment_call, first=TRUE) {
    gene_regex <- '(IG[HLK][VDJ]\\d+[-/\\w]*)'
    
    if (first) {
        r <- stringr::str_extract(segment_call, stringr::perl(gene_regex))
    } else {
        r <- stringr::str_extract_all(segment_call, stringr::perl(gene_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}


#' @rdname getAllele
#' @export
getFamily <- function(segment_call, first=TRUE) {
    family_regex <- '(IG[HLK][VDJ]\\d+)'
    
    if (first) {
        r <- stringr::str_extract(segment_call, stringr::perl(family_regex))
    } else {
        r <- stringr::str_extract_all(segment_call, stringr::perl(family_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}

