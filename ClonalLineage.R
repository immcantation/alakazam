#' Ig lineage reconstruction via maximum parsimony
#' 
#' @author     Jason Anthony Vander Heiden
#' @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @version    0.2.0
#' @date       2014.9.16

#### Imports ####
source_path <- dirname(sys.frame(1)$ofile)
source(file.path(source_path, "IOCore.R"), chdir=TRUE)
source(file.path(source_path, "SeqCore.R"), chdir=TRUE)


# Multifurcating tree
#ape::multi2di
# Read tree file into phylo object
#ape::read.tree
# Make/score trees with parsimony
#phangorn::parsimony


#' Load clones from clip tab delimited file
#'
#' @param   sample_name   sample file identifier
#' @param   sample_file   clip tab delimited file
#' @return  a data.frame of the processed clones
processClipClones <- function(sample_name, sample_file) {
    # Load clones
    cat('\n', sample_name, ' -> LOADING FILE ', format(Sys.time(), '%Y-%m-%d %H:%M'), '\n--------------------------------\n', sep='')
    clip_df <- loadClip(sample_file)
    bad_idx <- which(nchar(clip_df$SEQUENCE_GAP) != nchar(clip_df$GERMLINE_GAP_D_MASK))
    if (length(bad_idx) > 0) {
        cat('INFO> Removing', length(bad_idx), 
            'sequences with mismatched germline and sequence lengths.\n')
        clip_df <- clip_df[-bad_idx,]
    }
    
    clip_df$SAMPLE <- sample_name
    clip_df$SEQUENCE_GAP <- replaceSeqGaps(toupper(clip_df$SEQUENCE_GAP))
    clip_df$GERMLINE_GAP_D_MASK <- replaceSeqGaps(toupper(clip_df$GERMLINE_GAP_D_MASK))
    clip_df$CLONE <- paste(clip_df$SAMPLE, clip_df$CLONE, sep='-')
    
    len_sum <- ddply(clip_df, .(CLONE), summarize, 
                     LENS=length(unique(nchar(SEQUENCE_GAP))))
    bad_cln <- len_sum$CLONE[len_sum$LENS > 1]
    if (length(bad_cln) > 0) {
        cat('INFO> Removing', length(bad_cln), 
            'clones with mismatched sequence lengths.\n')
        clip_df <- subset(clip_df, !(CLONE %in% bad_cln))
    }
    
    # Remove duplicates and collapse fields
    cat(sample_name, '-> REMOVING DUPLICATES', format(Sys.time(), '%Y-%m-%d %H:%M'), '\n')
    clip_df$taxa <- clip_df$SEQUENCE_ID
    clip_df$seq <- clip_df$SEQUENCE_GAP
    dedup_df <- ddply(clip_df, .(CLONE),
                      restrictedCollapse, 
                      uniq_fields=UNIQ_FIELDS, 
                      text_fields=TEXT_FIELDS, 
                      num_fields=NUM_FIELDS)
    dedup_df <- subset(dedup_df, select=-c(taxa, seq))
    
    cat(sample_name, '-> COUNTING MUTATIONS', format(Sys.time(), '%Y-%m-%d %H:%M'), '\n')
    dedup_df$MUT_COUNT <- mapply(seq.distance, dedup_df$SEQUENCE_GAP, dedup_df$GERMLINE_GAP_D_MASK)
    dedup_df$MUT_FREQ <- dedup_df$MUT_COUNT / nchar(dedup_df$SEQUENCE_GAP)
    
    cat(sample_name, '-> DONE', format(Sys.time(), '%Y-%m-%d %H:%M'), '\n')
    return(list(clones=dedup_df, raw=subset(clip_df, select=-c(taxa,seq))))
}

processesSequences <- function() {
    # Replace gaps with Ns, removed ragged edges, and collapse duplicates
    clone_df <- replace.gaps(clone_df, outer.only=F)
    seq_masked <- maskSequenceEnds(clone_df$seq[clone_df$taxa != germline], 
                                   max_mask=max_mask, trim=F)
    clone_df$seq[clone_df$taxa != germline] <- seq_masked
    unique_df <- collapseDuplicates(clone_df, text_fields=text_fields, 
                                    num_fields=num_fields, exclude_taxa=germline)
}