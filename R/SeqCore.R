#' Common DNA, amino acid, and gene annotation operations for Alakazam
#' 
#' @author     Jason Anthony Vander Heiden
#' @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @version    0.2.0
#' @date       2014.9.18


#### Imports ####
#require(seqinr)
require(plyr)
require(stringr)

#### Constants ####

# Default colors
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


#### Nucleotide functions ####

#' Build nucleotide distance matrix
#'
#' @param     gap    the value to assign to gap (-, .) characters
#' @return    matrix of nucleotide character distance
getNucMatrix <- function(gap=Inf) {
    IUPAC_DNA <- list('A'='A', 'C'='C', 'G'='G', 'T'='T',
                      'M'=c('A','C'), 'R'=c('A','G'), 'W'=c('A','T'), 
                      'S'=c('C','G'), 'Y'=c('C','T'), 'K'=c('G','T'), 
                      'V'=c('A','C','G'), 'H'=c('A','C','T'), 
                      'D'=c('A','G','T'), 'B'=c('C','G','T'),
                      'N'=c('A','C','G','T'))
    
    sub.mat <- diag(18)
    colnames(sub.mat) <- rownames(sub.mat) <- c(names(IUPAC_DNA), c('-','.','?'))
    for (i in 1:length(IUPAC_DNA)) {
        for (j in i:length(IUPAC_DNA)) {
            sub.mat[i, j] <- sub.mat[j, i] <- any(IUPAC_DNA[[i]] %in% IUPAC_DNA[[j]])
        }
    }
    sub.mat[c('.','-'), c('.','-')] <- 1 
    sub.mat[c('.','-'), 1:15] <- 1 - gap 
    sub.mat[1:15, c('.','-')] <- 1 - gap
    
    return(1 - sub.mat)
}


#' Determine sequence distance
#'
#' @param     seq1       sequence as a string
#' @param     seq2       sequence as a string
#' @param     nuc_mat    nucleotide character distance matrix
#' @return    distance between seq1 and seq2
getSeqDistance <- function(seq1, seq2, nuc_mat=getNucMatrix(gap=Inf)) {
    # Convert string to character vector
    seq1 <- unlist(strsplit(seq1, ''))
    seq2 <- unlist(strsplit(seq2, ''))
    valid.idx <- (seq1 != '-' | seq2 != '-')
    seq1 <- seq1[valid.idx]
    seq2 <- seq2[valid.idx]
    # Calculate distance
    d <- sapply(1:length(seq1), function(x) { nuc_mat[seq1[x], seq2[x]] })
    indels <- sum(rle(d)$values == Inf)
    
    return(sum(d[is.finite(d)]) + indels)
}


#' Replace gap characters with Ns in nucleotide sequences
#'
#' @param     sequences     a character vector of sequences
#' @param     outer_only    if TRUE replace only leading and trailing gaps
#' @return    a modified sequence vector
maskSeqGaps <- function(sequences, outer_only=FALSE) {
    if (outer_only) {
        for (i in 1:length(sequences)) {
            head_match <- attr(regexpr('^[-\\.]+', sequences[i]), 'match.length')
            tail_match <- attr(regexpr('[-\\.]+$', sequences[i]), 'match.length')
            if (head_match > 0) { sequences[i] <- gsub('^[-\\.]+', paste(rep('N', head_match), collapse=''), sequences[i]) }
            if (tail_match > 0) { sequences[i] <- gsub('[-\\.]+$', paste(rep('N', tail_match), collapse=''), sequences[i]) }
        }
    } else {
        sequences <- gsub('[-\\.]', 'N', sequences)
    }
    
    return(sequences)
}


#' Replaces ragged leading and trailing edges of aligned sequences with Ns
#'
#' @param     sequences    a vector of sequence as strings
#' @param     max_mask     the maximum number of characters to mask
#'                         if NULL no threshold is set
#' @param     trim         if TRUE cut sequences rather than mask
#' @return    a modified sequence vector with masked or trimmed sequences
maskSeqEnds <- function(sequences, max_mask=NULL, trim=FALSE) {
    #cat('FUNCTION> maskSequenceEnds\n')
    
    # Find length of leading and trailing Ns
    left_lengths <- attr(regexpr('(^N*)', sequences, perl=T), 'capture.length')
    right_lengths <- attr(regexpr('(N*$)', sequences, perl=T), 'capture.length')
    
    # Mask to minimal inner sequence length
    left_mask <- min(max(left_lengths[, 1]), max_mask)
    right_mask <- min(max(right_lengths[, 1]), max_mask)
    seq_lengths <- nchar(sequences)
    if (trim) {
        sequences <- substr(sequences, left_mask + 1, seq_lengths - right_mask)
        #cat('INFO> Trimmed first', left_mask, 'and last', right_mask, 'characters\n')    
    } else {
        substr(sequences, 0, left_mask) <- paste(rep('N', left_mask), collapse='')
        substr(sequences, seq_lengths - right_mask + 1, seq_lengths + 1) <- paste(rep('N', right_mask), collapse='')
        #cat('INFO> Masked first', left_mask, 'and last', right_mask, 'characters\n')
    }
    
    return(sequences)
}


#' Collapse duplicate sequences and fields
#'
#' @param   data           a ChangeoClone object return by prepareClone
#' @param   text_fields    a vector of text columns to collapse 
#' @param   num_fields     a vector of numeric columns to collapse 
#' @param   nuc_mat        nucleotide character distance matrix
#' @return  modified ChangeoClone object with duplicate sequences removed and fields collapsed
collapseDuplicates <- function(data, text_fields=NULL, num_fields=NULL, 
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
                    dimnames=list(data[, "SEQUENCE_ID"], data[, "SEQUENCE_ID"]))
    for (i in 1:(nseq - 1)) {
        for (j in (i + 1):nseq) {
            d_mat[i, j] <- d_mat[j, i] <- getSeqDistance(data[i, "SEQUENCE_GAP"], 
                                                         data[j, "SEQUENCE_GAP"], 
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
    unique_df <- subset(data, SEQUENCE_ID %in% uniq_taxa)
    
    # Collapse duplicate sets and append entries to unique data.frame
    for (taxa in dup_taxa) {
        # Define row indices of identical sequences
        idx <- which(data[, "SEQUENCE_ID"] %in% taxa)
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
            seq_set <- unique(data[idx, "SEQUENCE_GAP"])
            unambig_len <- nchar(gsub("N", "", seq_set))
            tmp_df[, "SEQUENCE_GAP"] <- seq_set[which.max(unambig_len)]
        }
        
        # Add row to unique data.frame
        unique_df <- rbind(unique_df, tmp_df)
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

#' Get segment family calls
#'
#' @param     segment_call    a character vector containing segment calls delimited by commas
#' @param     first           if TRUE return only the first call
#'                            if FALSE return all calls delimited by commas
#' @return    a string containg the family calls
getFamily <- function(segment_call, first=TRUE) {
    family_regex <- '(IG[HLK][VDJ]\\d+)'
    
    if (first) {
        r <- str_extract(segment_call, perl(family_regex))
    } else {
        r <- str_extract_all(segment_call, perl(family_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}


#' Get segment gene calls
#'
#' @param     segment_call    a character vector containing segment calls delimited by commas
#' @param     first           if TRUE return only the first call
#'                            if FALSE return all calls delimited by commas
#' @return    a string containg the gene calls
getGene <- function(segment_call, first=TRUE) {
    gene_regex <- '(IG[HLK][VDJ]\\d+[-/\\w]*)'
    
    if (first) {
        r <- str_extract(segment_call, perl(gene_regex))
    } else {
        r <- str_extract_all(segment_call, perl(gene_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}


#' Get segment allele calls
#'
#' @param     segment_call    a string containing one or more segment calls delimited by commas
#' @param     first           if TRUE return only the first call
#'                            if FALSE return all calls delimited by commas
#' @return    a string containg the allele calls
getAllele <- function(segment_call, first=TRUE) {
    allele_regex <- '(IG[HLK][VDJ]\\d+[-/\\w]*[-\\*][\\.\\w]+)'
    
    if (first) {
        r <- str_extract(segment_call, perl(allele_regex))
    } else {
        r <- str_extract_all(segment_call, perl(allele_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}
