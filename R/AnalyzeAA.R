# Amino acid sequence properties

#' @include Alakazam.R
NULL


#### Constants ####

# Possible regions
VDJ_REGIONS <- c("CDR1", "CDR2", "CDR3", "FWR1", "FWR2", "FWR3")

# Hydropathy index scores
HYDROPATHY <- c(+1.8, -4.5, -3.5, -3.5, +2.5,
                -3.5, -3.5, -0.4, -3.2, +4.5,
                +3.8, -3.9, +1.9, +2.8, -1.6,
                -0.8, -0.7, -0.9, -1.3, +4.2)
names(HYDROPATHY) <- c("A", "R", "N", "D", "C",
                       "Q", "E", "G", "H", "I",
                       "L", "K", "M", "F", "P",
                       "S", "T", "W", "Y", "V")	


#' Translate nucleotide sequences to amino acids
#' 
#' \code{translateDNA} translates nucleotide sequences to AA using functions 
#' from seqinr.
#' 
#' @param   seq   DNA sequence (a string) to be converted to AAs.
#' @param   trim  boolean flag to remove 3 nts from both ends of seq
#'          (converts IMGT junction to CDR3 region).
#' 
#' @return  string, translated AA stretch
#' 
#' @seealso  \code{\link[seqinr]{translate}}.
#' 
#' @examples
#' library(alakazam)
#' # Load Change-O file
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' seq <- df$JUNCTION[1:3]
#' translateDNA(df$JUNCTION[1])
#' translateDNA(df$JUNCTION[1], trim=TRUE)
#' translateDNA("ACTGACTCGA")
#' 
#' @export
translateDNA <- function (seq, trim=FALSE) {
    # Function to translate a single string
    .translate <- function(x) {
        if (nchar(x) >= 3) {
            paste(seqinr::translate(unlist(strsplit(x, ""))), collapse="")
        } else {
            NA
        }
    }
    
    # Remove 3 nucleotides from each end
	# Eg,  "ACTGACTCGA" -> "GACT" (with "ACT" and "CGA" removed)
	if (trim) { seq <- substr(seq, 4, nchar(seq) - 3) }
    
    # Apply translation
    aa <- sapply(seq, .translate)

    return(aa)
}


#' Calculates GRAVY hydrophobicity score from a given amino acid sequence.
#'
#' @param    seq  string defining the amino acid sequence.
#' 
#' @return   GRAVY score for the sequence.
#'
#' @examples 
#' seq <- "CARDRSTPWRRGIASTTVRTSW"
#' gravy(seq)
#'
#' @export
gravy <- function(seq) {
    # Create character vector from string
    aa <- unlist(strsplit(as.character(seq), ""))
    # Ignore AAs that are "*" and "X"
    aa <- aa[aa != "X" & aa != "*"]
    return(sum(HYDROPATHY[aa])/length(aa))
}


# Count patterns
# 
# Counts the number of times a "pattern" occurs in "x", a string
#
# @param   x        a string (usually amino acids)
# @param   pattern  regular expression to be matched in string
#  
# @return  number of times the regular expression was found
countOccurrences <- function(x, pattern) {
    return(sum(gregexpr(pattern, x)[[1]] > 0))
}


#' Count amino acid patterns
#' 
#' Counts the fraction of times amino acid patterns occur in sequences.
#'
#' @param   region      character vector of either nucleotide or amino acid sequences.
#' @param   patterns    list of amino acid patterns to count in seq. If list is named,
#'                      names will become column names of output data.frame.
#' @param   nt		    boolean, TRUE if the sequences (or sequence) are DNA.
#' @param   trim        if \code{TRUE} remove the first and last codon/amino acids from each
#'                      sequence before calculating properties. If \code{FALSE} do
#'                      not modify input sequences.
#' @param   region_name   name of region to add as prefix to output column names.
#'  
#' @return  data.frame, fraction of times each amino acid pattern was found.
#' 
#' @examples 
#' region <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
#'             "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
#'             "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
#' patterns <- c("A","V","[LI]")
#' names(patterns) <- c("ARG","VAL","ISO_LEU")
#' countPatterns(region, patterns, nt=TRUE, trim=TRUE, region_name="CDR3")
#'             
#' @export
countPatterns <- function(region, patterns, nt=FALSE, trim=FALSE, region_name="REGION") {
    # Translate sequences if nucleotide
    if (nt) {
        # Check for sequences that are too short
        not_empty <- which(nchar(region) > 2)
        
        # Translate all regions from nt to aa
        region_aa <- rep("", length(region))
        region_aa[not_empty] <- sapply(region[not_empty], translateDNA, trim=trim)
    } else {
        region_aa <- region
    }
    # Calculate region lengths
    aa_length <- sapply(region_aa, nchar)
    # Count occurrence of each amino acid pattern for each sequence
    out_df <- data.frame(matrix(0, nrow=length(region_aa), ncol=length(patterns)))
    names(out_df) <- paste(region_name, names(patterns), sep="_")
    for(i in 1:length(patterns)) {
        out_df[, i] <- sapply(region_aa, countOccurrences, patterns[i]) / aa_length
    }
    return(out_df)
}


#' Calculates amino acid properties of regions
#'
#' Obtain data.frame of Region properties. Note this can be used for any region. 
#'
#' @param   region        character vector of either nucleotide or amino acid sequences.
#' @param   nt      	  boolean, TRUE if the sequences (or sequence) are DNA.
#' @param   trim          if \code{TRUE} remove the first and last codon/amino acids from each
#'                        sequence before calculating properties. If \code{FALSE} do
#'                        not modify input sequences.
#' @param   region_name   name of region to add as prefix to output column names.
#' 
#' @return  data.frame, with columns for each property
#' 
#' @examples 
#' region <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
#'             "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
#'             "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
#' regionProperties(region, nt=TRUE, trim=TRUE, region_name="CDR3")
#' 
#' @export
regionProperties <- function(region, nt=FALSE, trim=FALSE, region_name="REGION") {
    #nt=T
    #trim=T
    #region_name=NULL
    
    if (nt) {
        # Check for sequences that are too short
        not_empty <- which(nchar(region) > 2)
        #message(paste(length(not_empty), "Region sequences found."))
        
        # Translate all regions from nt to aa
        region_aa <- rep("", length(region))
        region_aa[not_empty] <- sapply(region[not_empty], translateDNA, trim=trim)
        #message("Region translated to amino acids.")
    } else {
        region_aa <- region
    }
    
    # Calculate region lengths
    aa_length <- sapply(region_aa, nchar)

    # GRAVY (Grand Average of Hydropathy) index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    # Produces NA if there is a stop codon.
    aa_gravy <- sapply(region_aa, gravy)
    
    # Count the fraction of aa that are positively charged
    aa_positive <- sapply(region_aa, countOccurrences, "[RK]") / aa_length
    
    # Count fraction of aa that are negatively charged
    aa_negative <- sapply(region_aa, countOccurrences, "[DE]") / aa_length
    
    # TODO: Why are we normalizing aliphatic index by length?
    # Aliphatic index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    n_ala <- sapply(region_aa, countOccurrences, "[A]")
    n_val <- sapply(region_aa, countOccurrences, "[V]")
    n_leu_ile <- sapply(region_aa, countOccurrences, "[LI]")
    a <- 2.9
    b <- 3.9
    aa_aliphatic = (n_ala + a*n_val + b*n_leu_ile) / aa_length

    # Count fraction of aa that are aromatic
    aa_aromatic <- sapply(region_aa, countOccurrences, "[FWHY]") / aa_length
    
    # Return the data.frame with amino acid properties
    out_df <- data.frame(AA_LENGTH=aa_length, GRAVY=aa_gravy, 
                         AA_POSITIVE=aa_positive, AA_NEGATIVE=aa_negative, 
                         ALIPHATIC=aa_aliphatic, AROMATIC=aa_aromatic)
    colnames(out_df) <- paste0(region_name, "_", colnames(out_df))
    
    return(out_df)
    
}