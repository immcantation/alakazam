# Amino acid sequence properties

#' @include Alakazam.R
NULL


#### Constants ####

# Hydropathy index scores
# Equivalent to Peptides::H$KyteDoolittle
HYDROPATHY <- c("A"=+1.8, 
                "R"=-4.5, 
                "N"=-3.5, 
                "D"=-3.5, 
                "C"=+2.5,
                "Q"=-3.5, 
                "E"=-3.5, 
                "G"=-0.4, 
                "H"=-3.2, 
                "I"=+4.5,
                "L"=+3.8, 
                "K"=-3.9, 
                "M"=+1.9, 
                "F"=+2.8, 
                "P"=-1.6,
                "S"=-0.8, 
                "T"=-0.7, 
                "W"=-0.9, 
                "Y"=-1.3, 
                "V"=+4.2)

# Bulkiness (side chain volume/length) index scores
# Zimmerman et al, 1968
# TODO:  I don't see this in the Peptides package. Maybe it's there, maybe it isn't?
BULKINESS <- c("A"=11.50, 
               "R"=14.28, 
               "N"=12.82, 
               "D"=11.68, 
               "C"=13.46,
               "Q"=14.45, 
               "E"=13.57, 
               "G"=3.40, 
               "H"=13.69, 
               "I"=21.40,
               "L"=21.40, 
               "K"=15.71, 
               "M"=16.25, 
               "F"=19.8, 
               "P"=17.43,
               "S"=9.47, 
               "T"=15.77, 
               "W"=21.67, 
               "Y"=18.03, 
               "V"=21.57)


#### General sequence functions ####

#' Build an asymetrical hydropathy distance matrix
#'
#' \code{getPropertyMatrix} returns an asymmetrical distance matrix for property
#' changes between amino acid characters.
#' 
#' @param    property   string defining the type of property to calculate the distance for.
#'                      One of \code{c("hydropathy", "bulkiness")}. \code{"hydropathy"}
#'                      returns a matrix of hydrophobicity differences according to the 
#'                      Kyte & Doolittle scale. \code{"bulkiness"} returns a distance 
#'                      matrix according to the Zimmerman et al bulkiness index.
#'                      
#' @return   A \code{matrix} of amino acid character distances with row and 
#'           column names indicating the character pair.
#' 
#' @seealso  Create a distance matrix for \code{\link{getSeqDistance}}.
#'           See \link{getAAMatrix} for an amino acid Hamming distance matrix.
#' 
#' @examples
#' getPropertyMatrix("hydro")
#' getPropertyMatrix("bulk")
#' 
#' @export
getPropertyMatrix <- function(property) {
    # Check arguments
    property <- match.arg(property, c("hydropathy", "bulkiness"))
    
    if (property == "hydropathy") {
        prop_data <- HYDROPATHY
    } else if (property == "bulkiness") {
        prop_data <- BULKINESS
    }
    
    n <- length(prop_data)
    # Define Hamming distance matrix
    dist_mat <- matrix(0, n + 2, n + 2)
    colnames(dist_mat) <- rownames(dist_mat) <- c(names(prop_data), "X", "-")
    for (i in 1:n) {
        for (j in i:n) {
            dist_mat[i, j] <- prop_data[j] - prop_data[i]
            dist_mat[j, i] <- prop_data[i] - prop_data[j]
        }
    }
    
    return(dist_mat)
}


#' Calculate distance between two sequences
#' 
#' \code{getPropertyDistance} calculates the distance (change) in chemical properties 
#' between two amino acid sequences, excluding ambiguous positions.
#'
#' @param    seq1       character string containing an amino acid or DNA sequence. As property 
#'                      changes may be asymmetrical, \code{seq1} defines the starting sequence.
#' @param    seq2       character string containing an amino acid or DNA sequence. As property 
#'                      changes may be asymmetrical, \code{seq2} defines the ending sequence.
#' @param    property   string defining the type of property to calculate the distance for.
#'                      One of \code{c("hydropathy", "bulkiness")}.
#'                      \itemize{
#'                        \item \code{"hydropathy"}  distance is hydrophobicity differences 
#'                                                   according to the Kyte & Doolittle scale. 
#'                        \item \code{"bulkiness"}   distance is side chain bulk differences 
#'                                                   according to the Zimmerman et al 
#'                                                   bulkiness index.
#'                      }
#' @param    nt         specify \code{TRUE} if the sequence(s) are DNA and require translation
#'                      before calculating distance.
#' @param    normalize  string defining the normalization method.  If \code{"length"} then the
#'                      distance is normalized by the number of shared unambiguous positions
#'                      (informative positions). If \code{"none"} then no normalization is 
#'                      performed.
#'
#' @return   Numerical distance for the given property from \code{seq1} to \code{seq2}.
#' 
#' @seealso  Raw distance matrix can get generated using \code{\link{getPropertyMatrix}}.
#'           
#' @examples
#' # Hydrophobicity of unambiguous DNA seqences
#' getPropertyDistance("ATGGCC", "ATGGGC", "hydro", nt=TRUE)
#' getPropertyDistance("ATGGCC", "ATGGGC", "hydro", nt=TRUE, normalize="none")
#' # Bulkiness of unambiguous DNA seqences
#' getPropertyDistance("ATGGCC", "ATGGGC", "bulk", nt=TRUE)
#' getPropertyDistance("ATGGCC", "ATGGGC", "bulk", nt=TRUE, normalize="none")
#' 
#' # Amino acid sequence with ambiguous (X) positions
#' getPropertyDistance("AYQXG", "ATXXG", "hydro")
#' getPropertyDistance("AYQXG", "ATXXG", "hydro", normalize="none")
#' 
#' 
#' @export
getPropertyDistance <- function(seq1, seq2, property, nt=FALSE, 
                                normalize=c("length", "none")) {
    # Check arguments
    property <- match.arg(property, c("hydropathy", "bulkiness"))
    normalize <- match.arg(normalize)
    
    # Check input
    if (nchar(seq1) != nchar(seq2)) {
        stop("seq1 and seq2 must be equal in length")
    }
    
    # Define property distances
    dist_mat <- getPropertyMatrix(property)

    # Translate
    if (nt) { 
        seq1 <- translateDNA(seq1) 
        seq2 <- translateDNA(seq2) 
    }
    
    # Convert string to character vector
    seq1 <- unlist(strsplit(seq1, ""), use.names=FALSE)
    seq2 <- unlist(strsplit(seq2, ""), use.names=FALSE)
    
    # Remove positions that are ambiguous in either sequence
    i <- !(seq1 %in% c("X", "-", ".")) & !(seq2 %in% c("X", "-", "."))
    seq1 <- seq1[i]
    seq2 <- seq2[i]
    n <- length(seq1)
    
    # Calculate distance
    d <- sapply(1:n, function(x) { dist_mat[seq1[x], seq2[x]] })
    d <- sum(d)
    
    if (normalize == "length") {
        d <- d / n
    }
    
    return(d)
}


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


#### Chemical property functions ####

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
    # If patterns are unnamed, make the names X1...Xn
    if(is.null(names(patterns))) { names(patterns) <- names(out_df) }
    # If region name, append to names of patterns
    if(region_name != '') { 
        names(out_df) <- paste(region_name, names(patterns), sep="_") 
    } else {
        names(out_df) <- names(patterns)
    }
    # Iterate over patterns
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