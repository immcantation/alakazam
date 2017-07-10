# Common DNA, amino acid, and gene annotation operations for Alakazam


#### Sequence distance functions ####

#' Build a DNA distance matrix
#'
#' \code{getDNAMatrix} returns a Hamming distance matrix for IUPAC ambiguous
#' DNA characters with modifications for gap, \code{c("-", ".")}, and missing, 
#' \code{c("?")}, character values.
#' 
#' @param    gap  value to assign to characters in the set \code{c("-", ".")}.
#' 
#' @return   A \code{matrix} of DNA character distances with row and column names 
#'           indicating the character pair. By default, distances will be either 0 
#'           (equivalent), 1 (non-equivalent or missing), or -1 (gap). 
#' 
#' @seealso  Creates DNA distance matrix for \link{seqDist}.
#'           See \link{getAAMatrix} for amino acid distances.
#' 
#' @examples
#' # Set gap characters to Inf distance
#' # Distinguishes gaps from Ns
#' getDNAMatrix()
#' 
#' # Set gap characters to 0 distance
#' # Makes gap characters equivalent to Ns
#' getDNAMatrix(gap=0)
#' 
#' @export
getDNAMatrix <- function(gap=-1) {
    # Define Hamming distance matrix
    sub_mat <- diag(18)
    colnames(sub_mat) <- rownames(sub_mat) <- c(names(IUPAC_DNA), c("-", ".", "?"))
    for (i in 1:length(IUPAC_DNA)) {
        for (j in i:length(IUPAC_DNA)) {
            sub_mat[i, j] <- sub_mat[j, i] <- any(IUPAC_DNA[[i]] %in% IUPAC_DNA[[j]])
        }
    }
    
    # Add gap characters
    sub_mat[c("-", "."), c("-", ".")] <- 1 
    sub_mat[c("-", "."), 1:15] <- 1 - gap 
    sub_mat[1:15, c("-", ".")] <- 1 - gap
    
    return(1 - sub_mat)
}

#' Build an AA distance matrix
#'
#' \code{getAAMatrix} returns a Hamming distance matrix for IUPAC ambiguous
#' amino acid characters.
#' 
#' @param    gap  value to assign to characters in the set \code{c("-", ".")}.
#' 
#' @return   A \code{matrix} of amino acid character distances with row and column names 
#'           indicating the character pair.
#' 
#' @seealso  Creates an amino acid distance matrix for \link{seqDist}.
#'           See \link{getDNAMatrix} for nucleotide distances.
#' 
#' @examples
#' getAAMatrix()
#' 
#' @export
getAAMatrix <- function(gap=0) {
  # Define Hamming distance matrix
  sub_mat <- diag(27)
  colnames(sub_mat) <- rownames(sub_mat) <- c(names(IUPAC_AA), c("*", "-", "."))
  for (i in 1:length(IUPAC_AA)) {
    for (j in i:length(IUPAC_AA)) {
      sub_mat[i, j] <- sub_mat[j, i] <- any(IUPAC_AA[[i]] %in% IUPAC_AA[[j]])
    }
  }
  
  # Add gap characters
  sub_mat[c("-", "."), c("-", ".")] <- 1 
  sub_mat[c("-", "."), c(1:27)] <- 1 - gap 
  sub_mat[c(1:27), c("-", ".")] <- 1 - gap
  
  return(1 - sub_mat)
}


#### Sequence manipulation functions ####

#' Translate nucleotide sequences to amino acids
#' 
#' \code{translateDNA} translates nucleotide sequences to amino acid sequences.
#' 
#' @param   seq     vector of strings defining DNA sequence(s) to be converted to translated.
#' @param   trim    boolean flag to remove 3 nts from both ends of seq
#'          (converts IMGT junction to CDR3 region).
#' 
#' @return  A vector of translated sequence strings.
#' 
#' @seealso  \code{\link[seqinr]{translate}}.
#' 
#' @examples
#' # Translate a single sequence
#' translateDNA("ACTGACTCGA")
#'
#' # Translate a vector of sequences
#' translateDNA(ExampleDb$JUNCTION[1:3])
#' 
#' # Remove the first and last codon from the translation
#' translateDNA(ExampleDb$JUNCTION[1:3], trim=TRUE)
#' 
#' @export
translateDNA <- function (seq, trim=FALSE) {
    # Function to translate a single string
    .translate <- function(x) {
        if (stri_length(x) >= 3) {
            stri_join(seqinr::translate(unlist(strsplit(x, "")), ambiguous=TRUE), 
                      collapse="")
        } else {
            NA
        }
    }
    
    # Remove 3 nucleotides from each end
    # Eg,  "ACTGACTCGA" -> "GACT" (with "ACT" and "CGA" removed)
    if (trim) { seq <- substr(seq, 4, stri_length(seq) - 3) }
    
    # Replace gaps with N
    seq <- gsub("[-.]", "N", seq)
    
    # Apply translation
    aa <- sapply(seq, .translate, USE.NAMES=FALSE)
    
    return(aa)
}


#' Masks gap characters in DNA sequences
#' 
#' \code{maskSeqGaps} substitutes gap characters, \code{c("-", ".")}, with \code{"N"} 
#' in a vector of DNA sequences.
#'
#' @param    seq         a character vector of DNA sequence strings.
#' @param    outer_only  if \code{TRUE} replace only contiguous leading and trailing gaps;
#'                       if \code{FALSE} replace all gap characters.
#'                       
#' @return   A modified \code{seq} vector with \code{"N"} in place of \code{c("-", ".")} 
#'           characters.
#' 
#' @seealso  See \link{maskSeqEnds} for masking ragged edges.
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


#' Masks ragged leading and trailing edges of aligned DNA sequences
#' 
#' \code{maskSeqEnds} takes a vector of DNA sequences, as character strings,
#' and replaces the leading and trailing characters with \code{"N"} characters to create 
#' a sequence vector with uniformly masked outer sequence segments.
#' 
#' @param    seq       a character vector of DNA sequence strings.
#' @param    max_mask  the maximum number of characters to mask. If set to 0 then
#'                     no masking will be performed. If set to \code{NULL} then the upper 
#'                     masking bound will be automatically determined from the maximum 
#'                     number of observed leading or trailing \code{"N"} characters amongst 
#'                     all strings in \code{seq}. 
#' @param    trim      if \code{TRUE} leading and trailing characters will be cut rather 
#'                     than masked with \code{"N"} characters.
#' @return   A modified \code{seq} vector with masked (or optionally trimmed) sequences.
#' 
#' @seealso   See \link{maskSeqGaps} for masking internal gaps.
#' 
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
    seq_lengths <- stri_length(seq)
    if (trim) {
        seq <- substr(seq, left_mask + 1, seq_lengths - right_mask)
    } else {
        substr(seq, 0, left_mask) <- paste(rep('N', left_mask), collapse='')
        substr(seq, seq_lengths - right_mask + 1, seq_lengths + 1) <- 
            paste(rep('N', right_mask), collapse='')
    }
    
    return(seq)
}


#' Remove duplicate DNA sequences and combine annotations
#' 
#' \code{collapseDuplicates} identifies duplicate DNA sequences, allowing for ambiguous 
#' characters, removes the duplicate entries, and combines any associated annotations.
#'
#' @param    data         data.frame containing Change-O columns. The data.frame 
#'                        must contain, at a minimum, a unique identifier column 
#'                        and a column containg a character vector of DNA sequences.
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing DNA sequences.
#' @param    text_fields  character vector of textual columns to collapse. The textual 
#'                        annotations of duplicate sequences will be merged into a single 
#'                        string with each unique value alphabetized and delimited by 
#'                        \code{sep}.
#' @param    num_fields   vector of numeric columns to collapse. The numeric annotations
#'                        of duplicate sequences will be summed. 
#' @param    seq_fields   vector of nucletoide sequence columns to collapse. The sequence 
#'                        with the fewest numer of non-informative characters will be 
#'                        retained. Where a non-informative character is one of 
#'                        \code{c("N", "-", ".", "?")}. Note, this is distinct from the 
#'                        \code{seq} parameter which is used to determine duplicates.
#' @param    add_count    if \code{TRUE} add the column \code{COLLAPSE_COUNT} that 
#'                        indicates the number of sequences that were collapsed to build 
#'                        each unique entry.
#' @param    ignore       vector of characters to ignore when testing for equality.
#' @param    sep          character to use for delimiting collapsed annotations in the 
#'                        \code{text_fields} columns. Defines both the input and output 
#'                        delimiter.
#' @param    verbose      if \code{TRUE} report the number input, discarded and output 
#'                        sequences; if \code{FALSE} process sequences silently.
#' @param    dry          dry run collapseDuplicates. Will only label the sequences.
#'                        
#' @return   A modified \code{data} data.frame with duplicate sequences removed and 
#'           annotation fields collapsed. If dry=TRUE, sequences will not be collapsed
#'           and the columns COLLAPSE_ID, COLLAPSE_CLASS and COLLAPSE_PASS will be
#'           added to the data.frame. Sequence with the same COLLAPSE_ID, belong to the 
#'           same collapsing group. COLLAPSE_CLASS will be unique, duplicated or ambiguous.
#'           COLLAPSE_PASS will be TRUE for the sequences that collapseDuplicates 
#'           (with dry=FALSE) would keep.
#'           
#' @details
#' \code{collapseDuplicates} identifies duplicate sequences in the \code{seq} column by
#' testing for character identity, with consideration of IUPAC ambiguous nucleotide codes. 
#' A cluster of sequences are considered duplicates if they are all equivalent, and no 
#' member of the cluster is equivalent to a sequence in a different cluster. 
#' 
#' Textual annotations, specified by \code{text_fields}, are collapsed by taking the unique
#' set of values within in each duplicate cluster and delimiting those values by \code{sep}.
#' Numeric annotations, specified by \code{num_fields}, are collapsed by summing all values 
#' in the duplicate cluster. Sequence annotations, specified by \code{seq_fields}, are 
#' collapsed by retaining the first sequence with the fewest number of N characters.
#' 
#' Columns that are not specified in either \code{text_fields}, \code{num_fields}, or 
#' \code{seq_fields} will be retained, but the value will be chosen from a random entry 
#' amongst all sequences in a cluster of duplicates.
#' 
#' An ambiguous sequence is one that can be assigned to two different clusters, wherein
#' the ambiguous sequence is equivalent to two sequences which are themselves 
#' non-equivalent. Ambiguous sequences arise due to ambiguous characters at positions that
#' vary across sequences, and are discarded along with their annotations. Thus, ambiguous
#' sequences are removed as duplicates of some sequence, but do not create a potential
#' false-positive annotation merger. Ambiguous sequences are not included in the 
#' \code{COLLAPSE_COUNT} annotation that is added when \code{add_count=TRUE}.
#' 
#' @seealso  Equality is tested with \link{seqEqual} and \link{pairwiseEqual}. 
#'           For IUPAC ambiguous character codes see \link{IUPAC_DNA}.
#'
#' @examples
#' # Example Change-O data.frame
#' db <- data.frame(SEQUENCE_ID=LETTERS[1:4],
#'                  SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
#'                  TYPE=c("IgM", "IgG", "IgG", "IgA"),
#'                  SAMPLE=c("S1", "S1", "S2", "S2"),
#'                  COUNT=1:4,
#'                  stringsAsFactors=FALSE)
#' 
#' # Annotations are not parsed if neither text_fields nor num_fields is specified
#' # The retained sequence annotations will be random
#' collapseDuplicates(db, verbose=TRUE)
#' 
#' # Unique text_fields annotations are combined into a single string with ","
#' # num_fields annotations are summed
#' # Ambiguous duplicates are discarded
#' collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
#'                    verbose=TRUE)
#'
#' # Use alternate delimiter for collapsing textual annotations
#' collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
#'                    sep="/", verbose=TRUE)
#' 
#' # Add count of duplicates
#' collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
#'                    add_count=TRUE, verbose=TRUE)
#' 
#' # Masking ragged ends may impact duplicate removal
#' db$SEQUENCE_IMGT <- maskSeqEnds(db$SEQUENCE_IMGT)
#' collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
#'                    add_count=TRUE, verbose=TRUE)
#'
#' @export
collapseDuplicates <- function(data, id="SEQUENCE_ID", seq="SEQUENCE_IMGT",
                               text_fields=NULL, num_fields=NULL, seq_fields=NULL,
                               add_count=FALSE, ignore=c("N", "-", ".", "?"), 
                               sep=",", verbose=FALSE, dry=FALSE) {
    # Verify column classes and exit if they are incorrect
    if (!is.null(text_fields)) {
        if (!all(sapply(subset(data, select=text_fields), is.character))) {
            stop("All text_fields columns must be of type 'character'")
        }
    }
    if (!is.null(num_fields)) {
        if (!all(sapply(subset(data, select=num_fields), is.numeric))) {
            stop("All num_fields columns must be of type 'numeric'")
        }
    }
    if (!is.null(seq_fields)) {
        if (!all(sapply(subset(data, select=seq_fields), is.character))) {
            stop("All seq_fields columns must be of type 'character'")
        }
    }
    seq_len <- stri_length(data[[seq]])
    if (any(seq_len != seq_len[1])) {
        warning("All sequences are not the same length")
    }
    
    # Define verbose reporting function
    .printVerbose <- function(n_total, n_unique, n_discard) {
        cat(" FUNCTION> collapseDuplicates\n", sep="")
        cat("    TOTAL> ", n_total, "\n", sep="")
        cat("   UNIQUE> ", n_unique, "\n", sep="")
        cat("COLLAPSED> ", n_total - n_unique - n_discard, "\n", sep="")
        cat("DISCARDED> ", n_discard, "\n", sep="")
        cat("\n")
    }
    
    # Define function to count informative positions in sequences
    .informativeLength <- function(x) {
        stri_length(gsub("[N\\-\\.\\?]", "", x, perl=TRUE))
    }
    
    # Initialize COLLAPSE_COUNT with 1 for each sequence
    if(add_count) {
        data[["COLLAPSE_COUNT"]] <- rep(1, nrow(data))
        num_fields <- c(num_fields, "COLLAPSE_COUNT")
    }
    
    # Initialize dry run columns 
    if (dry) {
        data$COLLAPSE_ID <- NA
        data$COLLAPSE_CLASS <- "duplicated"
        data$COLLAPSE_PASS <- TRUE
    }
    
    # Return input if there are no sequences to collapse
    nseq <- nrow(data)
    if (nseq <= 1) { 
        if (verbose) { .printVerbose(nseq, 1, 0) }
        return(data)
    }
    
    # Build distance matrix
    d_mat <- pairwiseEqual(data[[seq]])
    colnames(d_mat) <- rownames(d_mat) <- data[[id]]
    
    # Return input if no sequences are equal
    if (!any(d_mat[lower.tri(d_mat, diag=F)])) {
        if (verbose) { .printVerbose(nseq, nseq, 0) }
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
    
    if (dry & length(ambig_rows)>0) {
        data[["COLLAPSE_CLASS"]][ambig_rows] <- "ambiguous"
        data[["COLLAPSE_PASS"]][ambig_rows] <- FALSE
    }
    
    # Return single sequence if all sequence belong to ambiguous clusters
    if (discard_count == nrow(d_mat)) {
        inform_len <- .informativeLength(data[[seq]])
        selected <- which.max(inform_len)
        if (verbose) { .printVerbose(nseq, 0, discard_count - 1) }
        if (dry) {
            data[["COLLAPSE_PASS"]][-selected] <- FALSE
        } else {
            return(data[selected, ])
        }
    }
    
    # Exclude ambiguous sequences from clustering
    if (!dry & discard_count > 0) {
            d_mat <- d_mat[-ambig_rows, -ambig_rows]   
    }
    
    # Cluster remaining sequences into unique and duplicate sets
    dup_taxa <-  list()
    uniq_taxa <- character()
    done_taxa <- character()
    taxa_names <- rownames(d_mat)
    collapse_id <- 1
    
    for (taxa_i in 1:length(taxa_names)) {
        
        taxa <- taxa_names[taxa_i]
        
        # Skip taxa if previously assigned to a cluster
        if (taxa %in% done_taxa) { next }
        
        # Find all zero distance taxa
        idx <- which(d_mat[taxa, ])
        
        # Update vector of clustered taxa
        done_taxa <- c(done_taxa, taxa_names[idx])
        
        # Update collapse group
        if (dry) {
            data[["COLLAPSE_ID"]][idx] <- paste(data[["COLLAPSE_ID"]][idx], collapse_id, sep=",")
        }
        
        idx <- idx[idx %in% ambig_rows == FALSE]
        
        if (length(idx) == 1) {
            # Assign unique sequences to unique vector
            uniq_taxa <- append(uniq_taxa, taxa_names[idx])
            if (dry) {
                data[["COLLAPSE_CLASS"]][taxa_i] <- "unique"
                data[["COLLAPSE_PASS"]][taxa_i] <- TRUE
            }
        } else if (length(idx) > 1) {
            # Assign clusters of duplicates to duplicate list            
            dup_taxa <- c(dup_taxa, list(taxa_names[idx]))    
            if (dry) {
                # Keep COLLAPSE_PASS==TRUE for the sequence with the
                # larger number of informative positions 
                # (the first one if ties)
                data[["COLLAPSE_PASS"]][idx[-which.max(.informativeLength(data[[seq]][idx]))][1]] <- FALSE
            }
        } else {
            # Report error (should never occur)
            if (!dry) {
                stop("Error in distance matrix of collapseDuplicates")
            }
        }
        
        collapse_id <- collapse_id+1
    }
   
    if (dry) {
        data[["COLLAPSE_ID"]] <-  sub("^NA,","",data[["COLLAPSE_ID"]])
        return(data)
    }
    
    # Collapse duplicate sets and append entries to unique data.frame
    unique_list <- list(data[data[[id]] %in% uniq_taxa, ])
    for (taxa in dup_taxa) {
        # Define row indices of identical sequences
        idx <- which(data[[id]] %in% taxa)
        tmp_df <- data[idx[1], ]
        
        if (length(idx) > 1) {
            
            # Initialize with data from most informative sequence
            seq_set <- data[idx, c(id, seq)]
            inform_len <- .informativeLength(seq_set[[seq]])
            max_inform <- which.max(inform_len)[1] # if ties, pick first
            tmp_df <- data[idx[max_inform],]
            
            # Define set of text fields for row
            for (f in text_fields) {
                f_set <- na.omit(data[[f]][idx])
                if (length(f_set) > 0) {
                    f_set <- unlist(strsplit(f_set, sep))
                    f_set <- sort(unique(f_set))
                    f_val <- paste(f_set, collapse=sep)
                } else {
                    f_val <- NA
                }
                tmp_df[, f] <- f_val
            }
            
            # Sum numeric fields
            for (f in num_fields) {
                f_set <- na.omit(data[[f]][idx])
                if (length(f_set) > 0) { 
                    f_val <- sum(f_set) 
                } else { 
                    f_val <- NA 
                }
                tmp_df[, f] <- f_val
            }
            
            # Select sequence fields with fewest Ns
            for (f in seq_fields) {
                f_set <- na.omit(data[[f]][idx])
                if (length(f_set) > 0) {
                    f_len <- .informativeLength(f_set)
                    f_val <- f_set[which.max(f_len)]
                } else {
                    f_val <- NA
                }
                tmp_df[, f] <- f_val
            }
            
        }
        
        # Add row to unique list
        unique_list <- c(unique_list, list(tmp_df))
    }
    
    # Combine all rows into unique data.frame
    unique_df <- as.data.frame(bind_rows(unique_list))
    
    if (verbose) { .printVerbose(nseq, nrow(unique_df), discard_count) }
    return(unique_df)
}

#### Annotation functions ####

#' Extracts FWRs and CDRs from IMGT-gapped sequences
#' 
#' \code{extractVRegion} extracts the framework and complementarity determining regions of 
#' the V-segment for IMGT-gapped immunoglobulin (Ig) nucleotide sequences according to the 
#' IMGT numbering scheme.
#'
#' @param     sequences  character vector of IMGT-gapped nucleotide sequences.
#' @param     region     string defining the region(s) of the V-segment to extract. 
#'                       May be a single region or multiple regions (as a vector) from
#'                       \code{c("FWR1", "CDR1", "FWR2", "CDR2" ,"FWR3")}.  By default, all
#'                       regions will be returned.
#'                       
#' @return    If only one region is specified in the \code{region} argument, a character 
#'            vector of the extracted sub-sequences will be returned. If multiple regions 
#'            are specified, then a character matrix will be returned with columns 
#'            corresponding to the specified regions and a row for each entry in 
#'            \code{sequences}.
#' 
#' @seealso   IMGT-gapped region boundaries are defined in \link{IMGT_REGIONS}.
#' 
#' @references
#' \enumerate{
#'   \item  Lefranc M-P, et al. IMGT unique numbering for immunoglobulin and T cell 
#'            receptor variable domains and Ig superfamily V-like domains.
#'            Dev Comp Immunol. 2003 27(1):55-77.
#' }
#' 
#' @examples
#' # Assign example clone
#' clone <- subset(ExampleDb, CLONE == 3138)
#'
#' # Get all regions
#' extractVRegion(clone$SEQUENCE_IMGT)
#' 
#' # Get single region
#' extractVRegion(clone$SEQUENCE_IMGT, "FWR1")
#' 
#' # Get all CDRs
#' extractVRegion(clone$SEQUENCE_IMGT, c("CDR1", "CDR2"))
#' 
#' # Get all FWRs
#' extractVRegion(clone$SEQUENCE_IMGT, c("FWR1", "FWR2", "FWR3"))
#'
#' @export
extractVRegion <- function(sequences, region=c("FWR1", "CDR1", "FWR2", "CDR2" ,"FWR3")) {
    # Check region argument
    region <- match.arg(region, several.ok=TRUE)
    
    if (length(region) == 1) {
        sub_sequences <- substr(sequences, 
                                IMGT_REGIONS[[region]][1], 
                                IMGT_REGIONS[[region]][2])
    } else {
        sub_sequences <- sapply(region, function(x) substr(sequences, 
                                                           IMGT_REGIONS[[x]][1], 
                                                           IMGT_REGIONS[[x]][2]))
    }
    
    return(sub_sequences)
}


#### Rcpp distance wrappers ####

#' Calculate distance between two sequences
#' 
#' \code{seqDist} calculates the distance between two DNA sequences.
#'
#' @param    seq1      character string containing a DNA sequence.
#' @param    seq2      character string containing a DNA sequence.
#' @param    dist_mat  Character distance matrix. Defaults to a Hamming distance 
#'                     matrix returned by \link{getDNAMatrix}. If gap 
#'                     characters, \code{c("-", ".")}, are assigned a value of -1 
#'                     in \code{dist_mat} then contiguous gaps of any run length,
#'                     which are not present in both sequences, will be counted as a 
#'                     distance of 1. Meaning, indels of any length will increase
#'                     the sequence distance by 1. Gap values other than -1 will 
#'                     return a distance that does not consider indels as a special case.
#'
#' @return   Numerical distance between \code{seq1} and \code{seq2}.
#' 
#' @seealso  Nucleotide distance matrix may be built with 
#'           \link{getDNAMatrix}. Amino acid distance matrix may be built
#'           with \link{getAAMatrix}. Used by \link{pairwiseDist} for generating
#'           distance matrices. See \link{seqEqual} for testing sequence equivalence.
#'           
#' @examples
#' # Ungapped examples
#' seqDist("ATGGC", "ATGGG")
#' seqDist("ATGGC", "ATG??")
#' 
#' # Gaps will be treated as Ns with a gap=0 distance matrix
#' seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=0))
#' 
#' # Gaps will be treated as universally non-matching characters with gap=1
#' seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=1))
#' 
#' # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
#' seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))
#' 
#' # Gaps of equivalent run lengths are not counted as gaps
#' seqDist("ATG-C", "ATG-C", dist_mat=getDNAMatrix(gap=-1))
#'
#' # Overlapping runs of gap characters are counted as a single gap
#' seqDist("ATG-C", "AT--C", dist_mat=getDNAMatrix(gap=-1))
#' seqDist("A-GGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))
#' seqDist("AT--C", "AT--C", dist_mat=getDNAMatrix(gap=-1))
#' 
#' # Discontiguous runs of gap characters each count as separate gaps
#' seqDist("-TGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))
#' 
#' @export
seqDist <- function(seq1, seq2, dist_mat=getDNAMatrix()) {
    seqDistRcpp(seq1, seq2, dist_mat)
}


#' Calculate pairwise distances between sequences
#' 
#' \code{pairwiseDist} calculates all pairwise distance between a set of sequences.
#'
#' @param    seq       character vector containing a DNA sequences.
#' @param    dist_mat  Character distance matrix. Defaults to a Hamming distance 
#'                     matrix returned by \link{getDNAMatrix}. If gap 
#'                     characters, \code{c("-", ".")}, are assigned a value of -1 
#'                     in \code{dist_mat} then contiguous gaps of any run length,
#'                     which are not present in both sequences, will be counted as a 
#'                     distance of 1. Meaning, indels of any length will increase
#'                     the sequence distance by 1. Gap values other than -1 will 
#'                     return a distance that does not consider indels as a special case.
#'
#' @return   A matrix of numerical distance between each entry in \code{seq}. 
#'           If \code{seq} is a named vector, row and columns names will be added 
#'           accordingly.
#' 
#'           Amino acid distance matrix may be built with \link{getAAMatrix}. 
#'           Uses \link{seqDist} for calculating distances between pairs.
#'           See \link{pairwiseEqual} for generating an equivalence matrix.
#'           
#' @examples
#' # Gaps will be treated as Ns with a gap=0 distance matrix
#' pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
#'              dist_mat=getDNAMatrix(gap=0))
#' 
#' # Gaps will be treated as universally non-matching characters with gap=1
#' pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
#'              dist_mat=getDNAMatrix(gap=1))
#' 
#' # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
#' pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
#'              dist_mat=getDNAMatrix(gap=-1))
#' 
#' @export
pairwiseDist <- function(seq, dist_mat=getDNAMatrix()) {
    pairwiseDistRcpp(seq, dist_mat)
}