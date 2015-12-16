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
        scores <- pK$EMBOSS
        names(scores) <- rownames(scores)
    }
    
    return(scores)
}


#' Calculates GRAVY hydrophobicity score of amino acid sequences.
#'
#' @param    seq           vector of amino acid sequences.
#' @param    hydropathy    named vector defining hydropathy values for each amino acid.
#' 
#' @return   GRAVY score for the sequence.
#'
#' @examples 
#' seq <- "CARDRSTPWRRGIASTTVRTSW"
#' gravy(seq)
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
        # Ignore AAs that are "*" and "X"
        x <- x[!(x %in% c("X", "-", ".", "*"))]
        sum(hydropathy[x])/length(x)
    }
    aa_gravy <- sapply(aa, .gravy)
    return(aa_gravy)
}


#' Calculates aliphatic index of amino acid sequences.
#'
#' @param    seq vector of amino acid sequences.
#' 
#' @return   aliphatic index for the sequence.
#'
#' @examples 
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "ASTQMYVRT")
#' aliphatic(seq)
#'
#' @export
aliphatic <- function(seq) {
    n_ala <- countOccurrences(seq, "[A]")
    n_val <- countOccurrences(seq, "[V]")
    n_leu_ile <- countOccurrences(seq, "[LI]")
    aa_aliphatic = (n_ala + 2.9 * n_val + 3.9 * n_leu_ile) / nchar(seq)
    
    return(aa_aliphatic)
}


#' Calculates charge of amino acid sequences.
#'
#' @param    seq  vector of amino acid sequences.
#' @param    pH   environmental pH.
#' @param    pK   named vector defining pK values for each amino acid.
#' 
#' @return   charge of the sequence.
#'
#' @examples 
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "ASTQMYVRT") 
#' charge(seq)
#'
#' @export
charge <- function(seq, pH=7.4, pK=NULL) {
    # Get charge data
    if(is.null(pK)) {
        pK <- getPropertyData("charge")
    }
    
    # Calculate charge
    carg <- countOccurrences(seq, "R") * (1/(1 + 10^(1 * (pH - pK["R"]))))
    chis <- countOccurrences(seq, "H") * (1/(1 + 10^(1 * (pH - pK["H"]))))
    clys <- countOccurrences(seq, "K") * (1/(1 + 10^(1 * (pH - pK["K"]))))
    casp <- countOccurrences(seq, "D") * (-1/(1 + 10^(-1 * (pH - pK["D"]))))
    cglu <- countOccurrences(seq, "E") * (-1/(1 + 10^(-1 * (pH - pK["E"]))))
    ccys <- countOccurrences(seq, "C") * (-1/(1 + 10^(-1 * (pH - pK["C"]))))
    ctyr <- countOccurrences(seq, "Y") * (-1/(1 + 10^(-1 * (pH - pK["Y"]))))
    aa_charge <- carg + clys + chis + casp + cglu + ctyr + ccys
    
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
    return(sapply(gregexpr(pattern, x), function(y) { sum(y>0) }))
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
        region_aa[not_empty] <- translateDNA(region[not_empty], trim=trim)
    } else {
        region_aa <- region
    }
    # Calculate region lengths
    aa_length <- nchar(region_aa)
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
        out_df[, i] <- countOccurrences(region_aa, patterns[i]) / aa_length
    }
    return(out_df)
}


#' Calculates amino acid properties of a sequence region
#'
#' Obtain data.frame of region properties. Note this can be used for any region. 
#'
#' @param   data          \code{data.frame} containing sequence data.
#' @param   seq           \code{character} name of the column containing input 
#'                        sequences.
#' @param   nt      	  boolean, TRUE if the sequences (or sequence) are DNA.
#' @param   trim          if \code{TRUE} remove the first and last codon/amino acids from each
#'                        sequence before calculating properties. If \code{FALSE} do
#'                        not modify input sequences.
#' @param   label         name of region to add as prefix to output column names.
#' 
#' @return  data.frame, with columns for each property
#' 
#' @examples 
#' db <- data.frame("JUNCTION"=c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
#'                               "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
#'                               "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC"))
#' regionProperties(db, seq="JUNCTION", nt=TRUE, trim=TRUE, label="CDR3")
#' 
#' @export
regionProperties <- function(data, seq="JUNCTION", 
                             nt=FALSE, trim=FALSE, label=NULL) {
    
    if (length(seq)>1) {
        seq <- seq[1]
        warning("seq must be a vector of length 1. Only the first element will be used. ")
    }
    region <- as.vector(data[,seq])
    
    if (nt) {
        # Check for sequences that are too short
        not_empty <- which(nchar(region) > 2)
        #message(paste(length(not_empty), "Region sequences found."))
        
        # Translate all regions from nt to aa
        region_aa <- character(length(region))
        region_aa[not_empty] <- translateDNA(region[not_empty], trim=trim)
        #message("Region translated to amino acids.")
    } else {
        region_aa <- region
    }
    
    # Calculate region lengths
    aa_length <- nchar(region_aa)

    # GRAVY (Grand Average of Hydropathy) index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    # Produces NA if there is a stop codon.
    aa_gravy <- sapply(region_aa, gravy)
    
    # Count the fraction of aa that are positively charged
    aa_positive <- countOccurrences(region_aa, "[RK]") / aa_length
    
    # Count fraction of aa that are negatively charged
    aa_negative <- countOccurrences(region_aa, "[DE]") / aa_length
    
    # Aliphatic index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    aa_aliphatic = sapply(region_aa, aliphatic)

    # Count fraction of aa that are aromatic
    aa_aromatic <- countOccurrences(region_aa, "[FWHY]") / aa_length
    
    # Return the data.frame with amino acid properties
    out_df <- data.frame(AA_LENGTH=aa_length, GRAVY=aa_gravy, 
                         AA_POSITIVE=aa_positive, AA_NEGATIVE=aa_negative, 
                         ALIPHATIC=aa_aliphatic, AROMATIC=aa_aromatic)
    # If no label, use sequence column name
    if(is.null(label)) { label <- sequenceColumn }
    colnames(out_df) <- paste0(label, "_", colnames(out_df))
    
    return(cbind(db, out_df))
}