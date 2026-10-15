# The dense pairwiseEqual implementation of .collapseDuplicatesSingle, as it was
# before equalityGraph (fix_oom @ 5cd70e6). Kept only as a reference for
# testCollapseDuplicates.R, which checks that the sparse implementation returns
# identical results. Not used by the package.
.collapseDuplicatesDense <- function(data, id, seq, text_fields, num_fields, seq_fields,
                                      add_count, ignore, sep, dry, verbose) {
    # Define verbose reporting function
    .printVerbose <- function(n_total, n_unique, n_discard) {
        cat(" FUNCTION> collapseDuplicates\n", sep="")
        cat(" FIRST_ID> ", data[[id]][1], "\n", sep="")
        cat("    TOTAL> ", n_total, "\n", sep="")
        cat("   UNIQUE> ", n_unique, "\n", sep="")
        cat("COLLAPSED> ", n_total - n_unique - n_discard, "\n", sep="")
        cat("DISCARDED> ", n_discard, "\n", sep="")
        cat("\n")
    }
    
    # Define function to count informative positions in sequences.
    # Independent of ignore, which defines equality, not information content.
    .informativeLength <- function(x) {
        stri_length(gsub("[N\\-\\.\\?]", "", x, perl=TRUE))
    }
    
    # Initialize collapse_count with 1 for each sequence
    if(add_count) {
        data[["collapse_count"]] <- rep(1, nrow(data))
        num_fields <- c(num_fields, "collapse_count")
    }
    
    # Initialize dry run columns 
    if (dry) {
        data$collapse_id <- NA
        data$collapse_class <- NA
        data$collapse_pass <- TRUE
    }
    
    # Return input if there are no sequences to collapse
    nseq <- nrow(data)
    if (nseq <= 1) { 
        if (verbose) { .printVerbose(nseq, 1, 0) }
        if (dry) {
            data[['collapse_id']] <- 1
            data[['collapse_class']] <- "unique"
            data[['collapse_pass']] <- TRUE            
        }
        return(data)
    }
    
    # Build distance matrix
    exact_duplicates <- any(duplicated(data[[seq]]))
    d_mat <- pairwiseEqual(unique(data[[seq]]), ignore=ignore)
    colnames(d_mat) <- rownames(d_mat) <- unique(data[[seq]])
    n_uniqueseq <- nrow(d_mat)
    
    # Return input if no sequences are equal
    if (!any(rowSums(d_mat) > 1) & !exact_duplicates) {
        if (verbose) { .printVerbose(nseq, nseq, 0) }
        if (dry) {
            data[['collapse_id']] <- 1:nrow(data)
            data[['collapse_class']] <- "unique"
            data[['collapse_pass']] <- TRUE            
        }
        return(data)
    }        
    
    # Find sequences that will cluster ambiguously
    ambig_rows <- numeric()
    for (i in 1:n_uniqueseq) {
        idx <- which(d_mat[i, ])
        tmp_mat <- d_mat[idx, idx]
        if (!all(tmp_mat)) { 
            ambig_rows <- append(ambig_rows, i) 
        }
    }
    discard_count <- length(ambig_rows)
    
    # from ambiguous rows in d_mat to 
    # ambiguous rows in data
    data_ambig_rows <- data[[seq]] %in% rownames(d_mat)[ambig_rows]
    data_discard_count <- sum(data_ambig_rows)
    
    if (dry & length(ambig_rows)>0) {
        data[["collapse_class"]][data_ambig_rows] <- "ambiguous"
        data[["collapse_pass"]][data_ambig_rows] <- FALSE
    }
    
    # Return single sequence if all or all but one sequence belong to ambiguous clusters 
	if (nrow(data) - data_discard_count <= 1) {
        inform_len <- data.frame(list("inform_len"=.informativeLength(data[[seq]])))
        # For each ambiguous cluster, return the best sequence
        g <- igraph::simplify(igraph::graph_from_adjacency_matrix(d_mat))
        inform_len$clusters <- igraph::components(g)$membership[data[[seq]]]
        inform_len$select_id <- 1:nrow(inform_len)
        selected <- inform_len %>%
            dplyr::group_by(!!rlang::sym("clusters")) %>%
            dplyr::slice(which.max(!!rlang::sym("inform_len"))) %>%
            dplyr::ungroup() %>%
            dplyr::select(!!rlang::sym("select_id")) %>% unlist() 
        
        if (verbose) { .printVerbose(nseq, 0, discard_count - 1) }
        if (dry) {
            data[["collapse_id"]] <- inform_len$clusters
            data[["collapse_pass"]][selected] <- TRUE
        } else {
            return(data[selected, ])
        }
    }
    
    # Exclude ambiguous sequences from clustering
    if (!dry & discard_count > 0) {
            d_mat <- d_mat[-ambig_rows, -ambig_rows, drop = FALSE]  # 'drop = FALSE' to keep dataframe structure if only one sequence is left 
            data <- data[!data_ambig_rows,]
    }
    
    # Cluster remaining sequences into unique and duplicate sets
    dup_taxa <-  list()
    uniq_taxa <- character()
    done_taxa <- character()
    taxa_names <- rownames(d_mat)
    collapse_id <- 1
    
    for (taxa_i in 1:length(taxa_names)) {
        
        taxa <- taxa_names[taxa_i]
        data_taxa_i <- which(data[[seq]] %in% taxa)
        
        # Skip taxa if previously assigned to a cluster
        # or if ambiguous
        # (ambiguous taxa don't get their own collapse_id)
        if (taxa %in% done_taxa) { next }
        if (dry & taxa_i %in% ambig_rows) { next }
        
        # Find all zero distance taxa
        idx <- which(d_mat[taxa, ])
        # Translate from d_mat idx to data idx
        data_idx <- which(data[[seq]] %in% colnames(d_mat)[idx])
        
        # Update vector of clustered taxa
        done_taxa <- c(done_taxa, taxa_names[idx])
        
        # Update collapse group
        if (dry) {
            data[["collapse_id"]][data_idx] <- paste(data[["collapse_id"]][data_idx], collapse_id, sep=",")
        }
        
        if (dry) {
            #idx_copy <- idx
            data_idx_copy <- data_idx
            idx <- idx[idx %in% ambig_rows == FALSE]
            data_idx <- which(data[[seq]] %in% colnames(d_mat)[idx])
        }
        
        if (length(data_idx) == 1) {
            # Assign unique sequences to unique vector
            uniq_taxa <- append(uniq_taxa, taxa_names[idx])
            if (dry) {
                if (length(data_idx_copy)==1) {
                    ## 'truly' unique
                    data[["collapse_class"]][data_taxa_i] <- "unique"    
                } else {
                    ## unique after ambiguous removal
                    data[["collapse_class"]][data_taxa_i] <- "ambiguous_duplicate"
                }
                data[["collapse_pass"]][data_taxa_i] <- TRUE
            }
        } else if (length(data_idx) > 1) {
            # Assign clusters of duplicates to duplicate list            
            dup_taxa <- c(dup_taxa, list(taxa_names[idx]))    
            if (dry) {
                # Keep collpase_pass==TRUE for the sequence with the
                # larger number of informative positions 
                # (the first one if ties)
                max_info_idx <- which.max(.informativeLength(data[[seq]][data_idx]))[1]
                data[["collapse_class"]][data_idx] <- "duplicated"
                data[["collapse_pass"]][data_idx[-max_info_idx]] <- FALSE
            }
        } else {
            # Report error (should never occur)
                stop("Error in distance matrix of collapseDuplicates")
        }
        
        collapse_id <- collapse_id + 1
    }
   
    if (dry) {
        data[["collapse_id"]] <-  sub("^NA,","",data[["collapse_id"]])
        return(data)
    }
    
    # Collapse duplicate sets and append entries to unique data.frame
    unique_list <- list(data[data[[seq]] %in% uniq_taxa, ])
    for (taxa in dup_taxa) {
        # Define row indices of identical sequences
        idx <- which(data[[seq]] %in% taxa)
        tmp_df <- data[idx[1], ]

        if (length(idx) > 1) {

            # Initialize with data from most informative sequence
            seq_set <- data[idx, c(id, seq)]
            inform_len <- .informativeLength(seq_set[[seq]])
            max_inform <- which.max(inform_len)[1] # if ties, pick first
            tmp_df <- data[idx[max_inform], ]

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
