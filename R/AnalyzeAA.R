# Amino acid sequence properties

#' @include Alakazam.R
NULL

#### Chemical property functions ####

# Wrapper function to pull amino acid property vectors from seqinr
#
# @param    property  string defining the property. One of
#                     \code(c("hydropathy", "bulkiness", "polarity", "charge"))
#
# @return   A vector of scores for the given property with single
#           character amino acid labels.
#
# @examples
# getPropertyData("hydro")
# getPropertyData("bulk")
# getPropertyData("polar")
# getPropertyData("charge")
getPropertyData <- function(property){
    property <- match.arg(property, c("hydropathy", "bulkiness", "polarity" ,"charge"))
    
    if (property == "hydropathy") {
        # Kyte & Doolittle, 1982.
        data(aaindex, package="seqinr", envir=environment())
        scores <- aaindex[["KYTJ820101"]]$I
        names(scores) <- translateStrings(names(scores), AA_TRANS)
    } else if (property == "bulkiness") {
        # Zimmerman et al, 1968.
        data(aaindex, package="seqinr", envir=environment())
        scores <- aaindex[["ZIMJ680102"]]$I
        names(scores) <- translateStrings(names(scores), AA_TRANS)
    } else if (property == "polarity") {
        # Grantham, 1974
        data(aaindex, package="seqinr", envir=environment())
        scores <- aaindex[["GRAR740102"]]$I
        names(scores) <- translateStrings(names(scores), AA_TRANS)
    } else if (property == "charge") {
        # EMBOSS
        data(pK, package="seqinr", envir=environment())
        scores <- setNames(pK[["EMBOSS"]], rownames(pK))
    }
    
    return(scores)
}


#' Calculates the hydrophobicity of amino acid sequences
#'
#' \code{gravy} calculates the Grand Average of Hydrophobicity (GRAVY) index 
#' of amino acid sequences using the method of Kyte & Doolittle. Non-informative
#' positions are excluded, where non-informative is defined as any character in 
#' \code{c("X", "-", ".", "*")}.
#' 
#' @param    seq         vector of strings containing amino acid sequences.
#' @param    hydropathy  named numerical vector defining hydropathy index values for 
#'                       each amino acid, where names are single-letter amino acid 
#'                       character codes. If \code{NULL}, then the Kyte & Doolittle
#'                       scale is used.
#' 
#' @return   A vector of GRAVY scores for the sequence(s).
#' 
#' @references
#' \enumerate{
#'   \item  Kyte J, Doolittle RF. A simple method for displaying the hydropathic character 
#'            of a protein. J Mol Biol. 157, 105-32 (1982).
#' 
#' @seealso 
#' For additional hydrophobicity indices see \code{\link[seqinr]{aaindex}}.
#'
#' @examples
#' # Default scale
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
#' gravy(seq)
#'
#' # Use the Kidera et al, 1985 scores from the seqinr package
#' library(seqinr)
#' data(aaindex)
#' h <- aaindex[["KIDA850101"]]$I
#' # Rename the score vector to use single-letter codes
#' names(h) <- translateStrings(names(h), AA_TRANS)
#' # Calculate hydrophobicity
#' gravy(seq, hydropathy=h)
#' 
#' @export
gravy <- function(seq, hydropathy=NULL) {
    # Get hydrophobicity scores
    if (is.null(hydropathy)) {
        hydropathy <- getPropertyData("hydro")
    }
    # Create character vector from string
    aa <- strsplit(as.character(seq), "")
    
    # Function to translate a single string
    .gravy <- function(x) {
        # Ignore AAs that are non-informative
        x <- x[!(x %in% c("X", "-", ".", "*"))]
        sum(hydropathy[x]) / length(x)
    }
    
    aa_gravy <- sapply(aa, .gravy)
    
    return(aa_gravy)
}


#' Calculates the aliphatic index of amino acid sequences
#' 
#' \code{aliphatic} calculates the aliphatic index of amino acid sequences using 
#' the method of Ikai. Non-informative positions are excluded, where non-informative 
#' is defined as any character in \code{c("X", "-", ".", "*")}.
#'
#' @param    seq  vector of strings containing amino acid sequences.
#' 
#' @return   A vector of the aliphatic indices for the sequence(s).
#'
#' @references 
#' \enumerate{
#'   \item  Ikai AJ. Thermostability and aliphatic index of globular proteins. 
#'            J Biochem. 88, 1895-1898 (1980).
#' }
#' 
#' @examples 
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
#' aliphatic(seq)
#'
#' @export
aliphatic <- function(seq) {
    # Remove non-informative positions
    seq <- gsub("[X\\.\\*-]", "", seq)
    
    # Calculate aliphatic index
    ala <- countOccurrences(seq, "[A]")
    val <- countOccurrences(seq, "[V]")
    leu_ile <- countOccurrences(seq, "[LI]")
    aa_aliphatic = (ala + 2.9 * val + 3.9 * leu_ile) / nchar(seq, keepNA=TRUE)
    
    return(aa_aliphatic)
}


#' Calculates the net charge of amino acid sequences.
#'
#' \code{charge} calculates the net charge of amino acid sequences using 
#' the method of Moore, 1985, with exclusion of the C-terminus and N-terminus charges.
#' 
#' @param    seq        vector strings defining of amino acid sequences.
#' @param    pH         environmental pH.
#' @param    pK         named vector defining pK values for each charged amino acid,
#'                      where names are the single-letter amino acid character codes
#'                      \code{c("R", "H", "K", "D", "E", "C", "Y")}). If \code{NULL}, 
#'                      then the EMBOSS scale is used.
#' @param    normalize  if \code{TRUE} then divide the net charge of each amino acid 
#'                      sequence by the number of informative positions. Non-informative 
#'                      position are defined by the presence any character in 
#'                      \code{c("X", "-", ".", "*")}. If \code{FALSE} then return the raw
#'                      net charge.
#' 
#' @return   A vector of net charges for the sequence(s).
#'
#' @references
#' \enumerate{
#'   \item  Moore DS. Amino acid and peptide net charges: A simple calculational procedure. 
#'            Biochem Educ. 13, 10-11 (1985).
#'   \item  \url{http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html}
#'   }
#'   
#' @seealso 
#' For additional pK scales see \code{\link[seqinr]{pK}}.
#' 
#' @examples 
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT") 
#' # Normalized charge
#' charge(seq)
#' # Unnormalized charge
#' charge(seq, normalize=FALSE)
#' 
#' # Use the Murray et al, 2006 scores from the seqinr package
#' library(seqinr)
#' data(pK)
#' x <- setNames(pK[["Murray"]], rownames(pK))
#' # Calculate charge
#' charge(seq, pK=x, normalize=FALSE)
#'
#' @export
charge <- function(seq, pH=7.4, pK=NULL, normalize=TRUE) {
    # Get charge data
    if(is.null(pK)) {
        pK <- getPropertyData("charge")
    }
    
    # Calculate charge
    arg <- countOccurrences(seq, "R") * (1/(1 + 10^(1 * (pH - pK["R"]))))
    his <- countOccurrences(seq, "H") * (1/(1 + 10^(1 * (pH - pK["H"]))))
    lys <- countOccurrences(seq, "K") * (1/(1 + 10^(1 * (pH - pK["K"]))))
    asp <- countOccurrences(seq, "D") * (-1/(1 + 10^(-1 * (pH - pK["D"]))))
    glu <- countOccurrences(seq, "E") * (-1/(1 + 10^(-1 * (pH - pK["E"]))))
    cys <- countOccurrences(seq, "C") * (-1/(1 + 10^(-1 * (pH - pK["C"]))))
    tyr <- countOccurrences(seq, "Y") * (-1/(1 + 10^(-1 * (pH - pK["Y"]))))
    aa_charge <- arg + lys + his + asp + glu + tyr + cys
    
    if (normalize) {
        aa_charge <- aa_charge / nchar(gsub("[X\\.\\*-]", "", seq), keepNA=TRUE)
    }
    
    return(aa_charge)
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
    return(sapply(gregexpr(pattern, x), function(y) { sum(y > 0) }))
}


#' Count sequence patterns
#' 
#' \code{countPatterns} counts the fraction of times a set of character patterns occur 
#' in a set of sequences.
#'
#' @param   seq         character vector of either DNA or amino acid sequences.
#' @param   patterns    list of sequence patterns to count in each sequence. If the 
#'                      list is named, then names will be assigned as the column names of 
#'                      output data.frame.
#' @param   nt		    if \code{TRUE} then \code{seq} are DNA sequences and require 
#'                      translations before performing the pattern search.
#' @param   trim        if \code{TRUE} remove the first and last codon or amino acid from 
#'                      each sequence before the pattern search. If \code{FALSE} do
#'                      not modify the input sequences.
#' @param   label       string defining a label to add as a prefix to the output 
#'                      column names.
#'  
#' @return  A data.frame containing the fraction of times each sequence pattern was 
#'          found.
#' 
#' @examples 
#' seq <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
#'          "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
#'          "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
#' patterns <- c("A", "V", "[LI]")
#' names(patterns) <- c("ARG", "VAL", "ISO_LEU")
#' countPatterns(seq, patterns, nt=TRUE, trim=TRUE, label="CDR3")
#'             
#' @export
countPatterns <- function(seq, patterns, nt=FALSE, trim=FALSE, label="REGION") {
    # Translate sequences if nucleotide
    region_aa <- if (nt) { translateDNA(seq, trim=trim) } else { seq }
    
    # TODO: Check that NA is passed through correctly.
    # TODO: What is the proper length normalization? With or without non-informative position?
    
    # Calculate region lengths
    aa_length <- nchar(region_aa, keepNA=TRUE)
    # Count occurrence of each amino acid pattern for each sequence
    out_df <- data.frame(matrix(0, nrow=length(region_aa), ncol=length(patterns)))
    # If patterns are unnamed, make the names X1...Xn
    if(is.null(names(patterns))) { names(patterns) <- names(out_df) }
    # If region name, append to names of patterns
    if(label != '') { 
        names(out_df) <- paste(label, names(patterns), sep="_") 
    } else {
        names(out_df) <- names(patterns)
    }
    # Iterate over patterns
    for(i in 1:length(patterns)) {
        out_df[, i] <- countOccurrences(region_aa, patterns[i]) / aa_length
    }
    return(out_df)
}


#' Calculates amino acid properties for sequence data
#'
#' \code{regionProperties} calculates amino acid sequences properties, including
#' length, hydrophobicity, aliphatic index and net charge.
#'
#' @param   data          \code{data.frame} containing sequence data.
#' @param   seq           \code{character} name of the column containing input 
#'                        sequences.
#' @param   nt      	  boolean, TRUE if the sequences (or sequence) are DNA.
#' @param   trim          if \code{TRUE} remove the first and last codon/amino acids from each
#'                        sequence before calculating properties. If \code{FALSE} do
#'                        not modify input sequences.
#' @param   label         name of sequence region to add as prefix to output column names.
#' 
#' @return  A modified \code{data} data.frame with the following columns:
#'          \itemize{
#'            \item  \code{*_AA_LENGTH}:     number of amino acids.
#'            \item  \code{*_AA_GRAVY}:      grand average of hydrophobicity (GRAVY) index.
#'            \item  \code{*_AA_ALIPHATIC}:  aliphatic index.
#'            \item  \code{*_AA_CHARGE}:     normalized net charge.
#'          }
#'          
#'          Where \code{*} is the value from \code{label} or the name specified for 
#'          \code{seq} if \code{label=NULL}.
#' @details 
#' 
#' #' @references
#' \enumerate{
#'   \item  Kyte J, Doolittle RF. A simple method for displaying the hydropathic character 
#'            of a protein. J Mol Biol. 157, 105-32 (1982).
#'   \item  Ikai AJ. Thermostability and aliphatic index of globular proteins. 
#'            J Biochem. 88, 1895-1898 (1980).
#'   \item  Moore DS. Amino acid and peptide net charges: A simple calculational procedure. 
#'            Biochem Educ. 13, 10-11 (1985).
#' }
#' 
#' @seealso 
#' See \link{countPatterns} for counting the occurance of specific amino acid subsequences.
#' See \link{gravy}, \link{aliphatic} and \link{charge} for functions that calculate the
#' included properties individually.
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' df <- df[c(1,10,100), ]
#' 
#' prop <- regionProperties(df, seq="JUNCTION", nt=TRUE, trim=TRUE, label="CDR3")
#' prop[, c(1, 15:18)]
#' 
#' @export
regionProperties <- function(data, seq="JUNCTION", nt=FALSE, trim=FALSE, label=NULL) {
    # Check input
    if (length(seq) > 1) {
        stop("You may specify only one sequence column. seq must be a vector of length 1.")
    }
    check <- checkColumns(data, seq)
    if (check != TRUE) { stop(check) }
    
    # Get sequence vector and translate if required
    region <- as.character(data[[seq]])
    region_aa <- if (nt) { translateDNA(region, trim=trim) } else { region }

    # Calculate region lengths
    aa_length <- nchar(region_aa, keepNA=TRUE)

    # Hydrophobicity
    aa_gravy <- gravy(region_aa)

    # Aliphatic index
    aa_aliphatic <- aliphatic(region_aa)
    
    # Net charge
    aa_charge <- charge(region_aa)
    
    # TODO: move these into a separate amino acid composition function
    # TODO: Should only include informative positions in the denominator
    # Count the fraction of aa that are positively charged
    #aa_positive <- countOccurrences(region_aa, "[RK]") / aa_length
    # Count fraction of aa that are negatively charged
    #aa_negative <- countOccurrences(region_aa, "[DE]") / aa_length
    # Count fraction of aa that are aromatic
    #aa_aromatic <- countOccurrences(region_aa, "[FWHY]") / aa_length
    
    # Return the data.frame with amino acid properties
    out_df <- data.frame(AA_LENGTH=aa_length, 
                         AA_GRAVY=aa_gravy, 
                         AA_ALIPHATIC=aa_aliphatic,
                         AA_CHARGE=aa_charge)
    
    # If no label, use sequence column name
    if (is.null(label)) { label <- seq }
    colnames(out_df) <- paste0(label, "_", colnames(out_df))
    
    return(cbind(data, out_df))
}