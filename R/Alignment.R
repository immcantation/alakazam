#' Calculate junction region alignment properties
#'
#' \code{junctionAlignment} determines the number of deleted germline nucleotides in the 
#' junction region and the number of V gene and J gene nucleotides in the CDR3.
#'
#' @param   data                \code{data.frame} containing sequence data.
#' @param   germline_db         reference germline database for the V, D and J genes.
#'                              in \code{data}
#' @param   v_call              V gene assignment column.
#' @param   d_call              D gene assignment column.
#' @param   j_call              J gene assignment column.
#' @param   v_germline_start    column containing the start position of the alignment 
#'                              in the V reference germline.
#' @param   v_germline_end      column containing the end position of the alignment in the 
#'                              V reference germline.
#' @param   d_germline_start    column containing the start position of the alignment 
#'                              in the D reference germline.
#' @param   d_germline_end      column containing the start position of the alignment 
#'                              in the D reference germline.
#' @param   j_germline_start    column containing the start position of the alignment 
#'                              in the J reference germline.
#' @param   j_germline_end      column containing the start position of the alignment 
#'                              in the J reference germline.
#' @param   np1_length          combined length of the N and P regions between the 
#'                              V and D regions (heavy chain) or V and J regions (light chain).      
#' @param   np2_length          combined length of the N and P regions between the 
#'                              D and J regions (heavy chain).            
#' @param   junction            column containing the junction sequence.
#' @param   junction_length     column containing the length of the junction region in nucleotides.
#' @param   sequence_alignment  column containing the aligned sequence.
#' 
#' @return  A modified input \code{data.frame} with the following additional columns storing 
#'          junction alignment information:
#'          \enumerate{
#'              \item  \code{v_germline_deleted_3}:  number of 3' V germline nucleotides deleted.
#'              \item  \code{d_germline_deleted_5}:  number of 5' D germline nucleotides deleted.
#'              \item  \code{d_germline_deleted_3}:  number of 3' D germline nucleotides deleted.
#'              \item  \code{j_germline_deleted_5}:  number of 5' J germline nucleotides deleted.
#'              \item  \code{v_cdr3_length}:         number of sequence_alignment V nucleotides in the CDR3.
#'              \item  \code{j_cdr3_length}:         number of sequence_alignment J nucleotides in the CDR3.
#'          }
#' 
#' @examples
#' germline_db <- list(
#' "IGHV3-11*05"="CAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTCAAGCCTGGAGGGTCCCTGAGACT
#' CTCCTGTGCAGCCTCTGGATTCACCTTC............AGTGACTACTACATGAGCTGGATCCGCCAGGCTCCAG
#' GGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTTACACAAACTACGCAGACTCTGTGAAG
#' ...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGA
#' CACGGCCGTGTATTACTGTGCGAGAGA",
#' "IGHD3-10*01"="GTATTACTATGGTTCGGGGAGTTATTATAAC",
#' "IGHJ5*02"="ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
#' )
#' 
#' db <- junctionAlignment(SingleDb, germline_db)
#'
#' @export
junctionAlignment <- function(data, germline_db, 
                                  v_call="v_call",
                                  d_call="d_call",
                                  j_call="j_call",
                                  v_germline_start="v_germline_start",
                                  v_germline_end="v_germline_end",
                                  d_germline_start="d_germline_start",
                                  d_germline_end="d_germline_end",
                                  j_germline_start="j_germline_start",
                                  j_germline_end="j_germline_end",
                                  np1_length="np1_length",
                                  np2_length="np2_length",
                                  junction="junction",
                                  junction_length="junction_length",
                                  sequence_alignment="sequence_alignment") {
    
    # Check input
    check <- checkColumns(data, 
                          c(v_call, d_call, j_call, 
                          v_germline_start, v_germline_end,
                          d_germline_start, d_germline_end,
                          j_germline_start, j_germline_end,
                          np1_length, np2_length,
                          junction, junction_length,
                          sequence_alignment))
    if (check != TRUE) { stop(check) }
    
    # Get deletions
    for (i in 1:nrow(data))  {
        v_dels <- countDeleted(data[i,], 
                               allele_call=v_call, germline_start=v_germline_start, germline_end=v_germline_end, 
                               germline_db=germline_db, junction=junction, junction_length=junction_length, 
                               sequence_alignment=sequence_alignment)
        d_dels <- countDeleted(data[i,], 
                               allele_call=d_call, germline_start=d_germline_start, germline_end=d_germline_end, 
                               germline_db=germline_db, junction=junction, junction_length=junction_length, 
                                sequence_alignment=sequence_alignment)
        j_dels <- countDeleted(data[i,], 
                               allele_call=j_call, germline_start=j_germline_start, germline_end=j_germline_end, 
                               germline_db=germline_db, junction=junction, junction_length=junction_length, 
                               sequence_alignment=sequence_alignment)
        data[['v_germline_deleted_3']][i] <- v_dels[2]
        data[['d_germline_deleted_5']][i] <- d_dels[1]
        data[['d_germline_deleted_3']][i] <- d_dels[2]
        data[['j_germline_deleted_5']][i] <- j_dels[1]
        data[['v_cdr3_length']][i] <- v_dels[3]
        data[['j_cdr3_length']][i] <- j_dels[3]
    }
    
    return(data)
}

# Junction alignment helper
#
# Report the number of deleted germline nucleotides in the alignment
#
# @param    db_row               one row from a Rearrangement database.
# @param    allele_call          column containing gene assignments.
# @param    germline_start       column containing the start position of the alignment in the reference germline.
# @param    germline_end         column containing the end position of the alignment in the reference germline.
# @param    germline_db          reference germline database for the V, D and J genes.
# @param    junction             column containing the junction sequence.
# @param    junction_length      column containing the length of the  junction region in nucleotides.
# @param    sequence_alignment   column containing the aligned sequence.
# 
# @return   Alignment deletions
countDeleted <- function(db_row, allele_call, germline_start, germline_end, 
                         germline_db, junction, junction_length,
                         sequence_alignment) {
    # db_row: one row from data
    # allele_call: one of v,d,j
    # germline_db: the reference germline database used to assign genes. 
    allele <- getAllele(db_row[[allele_call]], first=T)
    deleted <- c(NA, NA, NA)
    
    # Check for valid allele information
    if (is.na(allele)) { 
        return(deleted) 
    }
    # Check for allele in reference germlines
    tryCatch(germline <- germline_db[[allele]],
             error=function(e) { stop(allele, " not found in germline_db.") })
    
    allele_germline_start <- as.numeric(db_row[[germline_start]])
    allele_germline_end <- as.numeric(db_row[[germline_end]])
    
    germline_head <- stringi::stri_sub(germline, 1, allele_germline_start - 1)
    deleted_head <- nchar(gsub("\\.", "", germline_head))
    
    germline_tail <- stringi::stri_sub(germline, allele_germline_end+1, nchar(germline))
    deleted_tail <- nchar(gsub("\\.", "", germline_tail))
    
    deleted[1] <- deleted_head
    deleted[2] <- deleted_tail
    
    if (is.na(db_row[[junction]])) {
        warning("NA junction found.")
        return (deleted)
    }
    if (!db_row[[junction_length]]>6) {
        message("Junction length <= 6.")
        return (deleted)
    }
    
    junction_len <- db_row[[junction_length]]
    junction_start <- 310
    # junction_end <- junction_start + junction_len - 1
    
    # get aligned junction end (counting gaps)
    seq_aln <- s2c(db_row[[sequence_alignment]]) != "-"
    seq_aln[1:junction_start-1] <- 0
    junction_end <- which(cumsum(seq_aln[1:length(seq_aln)]) > junction_len)[1] - 1
    
    # For V and J alleles, calculate number of nt in the CDR3
    germ_cdr3_length <- NA
    if (grepl("[Vv]", allele)) {
        last_cdr3_pre_np <- db_row[[germline_end]] - db_row[[germline_start]] + 1 
        first_cdr3_pre_np <- junction_start + 3   # without conserved 
        # len <- last_cdr3_pre_np - first_cdr3_pre_np + 1
        #germ_seq <- stringi::stri_sub(germline, db_row[[germline_end]]+1-len, db_row[[germline_end]] )
        germ_seq <- stringi::stri_sub(db_row[[sequence_alignment]], first_cdr3_pre_np, last_cdr3_pre_np )
        germ_cdr3_length <- nchar(gsub("[\\.-]", "", germ_seq))
    } else if (grepl("[Jj]", allele))  {
        j_aln_len <- db_row[[germline_end]] - db_row[[germline_start]] + 1 
        # germ_seq <- stringi::stri_sub(germline, db_row[[germline_start]], db_row[[germline_end]]-j_tail)
        germ_seq <- stringi::stri_sub(db_row[[sequence_alignment]], 
                                      nchar(db_row[[sequence_alignment]]) - j_aln_len + 1,
                                      junction_end - 3)
        germ_cdr3_length <- nchar(gsub("-", "", germ_seq))
    } 
    
    deleted <- c(deleted_head, deleted_tail, germ_cdr3_length)
    return(deleted)
}
