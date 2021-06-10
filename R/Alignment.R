# Alignment properties

#' @param    data          Rearrangement database
#' @param    germline_db Reference germline database, for the V,D and J alleles
#'                       in \code{data}
#' @param    v_call      Name of the column containing V-segment allele assignments.
#' @param    d_call      Name of the column containing D-segment allele assignments.
#' @param    j_call      Name of the column containing J-segment allele assignments.
#' @param    v_germline_start  Name of the column containing the number of the starting
#'                      position of the alignment in the V reference germline 
#'                      in \code{germline_db}.
#' @param    d_germline_start  Name of the column containing the number of the starting
#'                      position of the alignment in the D reference germline 
#'                      in \code{germline_db}.
#' @param    d_germline_start  Name of the column containing the number of the starting
#'                      position of the alignment in the J reference germline 
#'                      in \code{germline_db}.
#' @param    junction  Name of the column containing the junction sequence.
#' @param    sequence_alignment  Name of the column containing the aligned sequence.
#' @return   Six new columns are added to \code{data}:
#' \enumerate{
#'   \item v_germline_deleted_3 Number of 3' V germline nucleotides deleted
#'   \item d_germline_deleted_5 Number of 5' D germline nucleotides deleted
#'   \item d_germline_deleted_3 Number of 3' D germline nucleotides deleted
#'   \item j_germline_deleted_5 Number of 5' J germline nucleotides deleted
#'   \item v_cdr3_length Number of sequence_alignment V nucleotides in the CDR3
#'   \item j_cdr3_length Number of sequence_alignment J nucleotides in the CDR3
#' }
#' 
#' @examples
#' data(oneseq_db)
#' germline_db <- list(
#' "IGHV3-11*05"="CAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTCAAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTC............AGTGACTACTACATGAGCTGGATCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGT......AGTAGTTACACAAACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCGAGAGA",
#' "IGHD3-10*01"="GTATTACTATGGTTCGGGGAGTTATTATAAC",
#' "IGHJ5*02"="ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
#' )
#' oneseq_db <- calcJunctionAlignment(oneseq_db, germline_db)
#' oneseq_db %>% select(contains("deleted"))
#' @export
calcJunctionAlignment <- function(data, germline_db, 
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
    
    check <- checkColumns(data, 
                          c(v_call, d_call, j_call, 
                          v_germline_start, v_germline_end,
                          d_germline_start, d_germline_end,
                          j_germline_start, j_germline_end,
                          np1_length, np2_length,
                          junction, junction_length,
                          sequence_alignment))
    if (check != TRUE) { stop(check) }
    
    # db_row: one row from data
    # allele_call: one of v,d,j
    # germline_db: the reference germline database used to assign genes. 
    .countDeleted <- function(db_row, allele_call, germline_start,germline_end, 
                              germline_db, junction, junction_length,
                              sequence_alignment,
                              v_germline_start, 
                              v_germline_end) {
        allele <- getAllele(db_row[[allele_call]], first=T)
        deleted <- c(NA,NA,NA)

        if (!is.na(allele)) {
            
            tryCatch(germline <- germline_db[[allele]],
                     error=function(e) {
                         stop(allele ," not found in germline_db.")
                     }
            )
            allele_germline_start <- as.numeric(db_row[[germline_start]])
            allele_germline_end <- as.numeric(db_row[[germline_end]])
            
            germline_head <- stringi::stri_sub(germline,1, allele_germline_start-1)
            deleted_head <- nchar(gsub("\\.","",germline_head))
            
            germline_tail <- stringi::stri_sub(germline,allele_germline_end+1,nchar(germline))
            deleted_tail <- nchar(gsub("\\.","",germline_tail))
            
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
            junction_end <- which(cumsum(seq_aln[1:length(seq_aln)])>junction_len)[1] - 1
            
            # For V and J alleles, calculate number of nt in the CDR3
            germ_cdr3_length <- NA
            if (grepl("[Vv]", allele)) {
                last_cdr3_pre_np <- db_row[[germline_end]] - db_row[[germline_start]] + 1 
                first_cdr3_pre_np <- junction_start + 3   # without conserved 
                len <- last_cdr3_pre_np - first_cdr3_pre_np + 1
                #germ_seq <- stringi::stri_sub(germline, db_row[[germline_end]]+1-len, db_row[[germline_end]] )
                germ_seq <- stringi::stri_sub(db_row[[sequence_alignment]], first_cdr3_pre_np, last_cdr3_pre_np )
                germ_cdr3_length <- nchar(gsub("[\\.-]","",germ_seq))
            } else if (grepl("[Jj]", allele))  {
                j_aln_len <- db_row[[germline_end]] - db_row[[germline_start]] + 1 
                # germ_seq <- stringi::stri_sub(germline, db_row[[germline_start]], db_row[[germline_end]]-j_tail)
                germ_seq <- stringi::stri_sub(db_row[[sequence_alignment]], 
                                              nchar(db_row[[sequence_alignment]]) - j_aln_len + 1,
                                              junction_end - 3)
                germ_cdr3_length <- nchar(gsub("-","",germ_seq))
            } 
            
            deleted <- c(deleted_head, deleted_tail, germ_cdr3_length)
        } 
        deleted
    }
    
    for (i in 1:nrow(data))  {
        v_dels <- .countDeleted(data[i,], v_call, v_germline_start, v_germline_end, germline_db, junction, junction_length, sequence_alignment, v_germline_start, v_germline_end)
        d_dels <- .countDeleted(data[i,], d_call, d_germline_start, d_germline_end, germline_db, junction, junction_length, sequence_alignment, v_germline_start, v_germline_end)
        j_dels <- .countDeleted(data[i,], j_call, j_germline_start, j_germline_end, germline_db, junction,  junction_length, sequence_alignment, v_germline_start, v_germline_end)
        data[['v_germline_deleted_3']][i] <- v_dels[2]
        data[['d_germline_deleted_5']][i] <- d_dels[1]
        data[['d_germline_deleted_3']][i] <- d_dels[2]
        data[['j_germline_deleted_5']][i] <- j_dels[1]
        data[['v_cdr3_length']][i] <- v_dels[3]
        data[['j_cdr3_length']][i] <- j_dels[3]
        
    }
    data
}
