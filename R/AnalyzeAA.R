# Amino acid sequence properties

#' @include Alakazam.R
NULL


#### Constants ####

# Amino acid abbreviation translation table
AA_TRANS <- c("A"="Ala",
              "R"="Arg",
              "N"="Asn",
              "D"="Asp",
              "C"="Cys",
              "Q"="Gln",
              "E"="Glu",
              "G"="Gly",
              "H"="His",
              "I"="Ile",
              "L"="Leu",
              "K"="Lys",
              "M"="Met",
              "F"="Phe",
              "P"="Pro",
              "S"="Ser",
              "T"="Thr",
              "W"="Trp",
              "Y"="Tyr",
              "V"="Val")

# Hydropathy index scores
# TODO: seqinr::aaindex[[151]]$I

# Bulkiness (side chain volume/length) index scores
# Zimmerman et al, 1968
# TODO: seqinr::aaindex[[399]]$I

# Polarity index scores
# Grantham, 1974
# TODO: seqinr::aaindex[[111]]$I


#### General sequence functions ####

#' Build an asymetrical hydropathy distance matrix
#'
#' \code{getPropertyMatrix} returns an asymmetrical distance matrix for property
#' changes between amino acid characters.
#' 
#' @param    property   string defining the type of property to calculate the distance for.
#'                      One of \code{c("hydropathy", "bulkiness", "polarity")}. 
#'                      \itemize{
#'                        \item \code{"hydropathy"} returns a matrix of hydrophobicity 
#'                                                  differences according to the 
#'                                                  Kyte & Doolittle scale. 
#'                        \item \code{"bulkiness"}  side chain bulkiness to the 
#'                                                  Zimmerman et al bulkiness index.
#'                        \item \code{"polarity"}   side chain polarity according to the
#'                                                  Grantham scale.
#'                      }
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
    property <- match.arg(property, c("hydropathy", "bulkiness", "polarity"))
    
    if (property == "hydropathy") {
        prop_data <- HYDROPATHY
    } else if (property == "bulkiness") {
        prop_data <- BULKINESS
    } else if (property == "polarity") {
        prop_data <- POLARITY
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


#' Calculate chemical property distance between two sequences
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
#'                        \item \code{"polarity"}    distince is side chain polarity differences
#'                                                   according to the Grantham scale.
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
#' getPropertyDistance("ATGGCCCGA", "ATGGGCCGT", "hydro", nt=TRUE)
#' getPropertyDistance("ATGGCCCGA", "ATGGGCCGT", "hydro", nt=TRUE, normalize="none")
#' # Bulkiness of unambiguous DNA seqences
#' getPropertyDistance("ATGGCCCGA", "ATGGGCCGT", "bulk", nt=TRUE)
#' getPropertyDistance("ATGGCCCGA", "ATGGGCCGT", "bulk", nt=TRUE, normalize="none")
#' # Polarity of unambiguous DNA seqences
#' getPropertyDistance("ATGGCCCGA", "ATGGGCCGT", "polar", nt=TRUE)
#' getPropertyDistance("ATGGCCCGA", "ATGGGCCGT", "polar", nt=TRUE, normalize="none")
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
    property <- match.arg(property, c("hydropathy", "bulkiness", "polarity"))
    normalize <- match.arg(normalize)
    
    # Check input
    if (nchar(seq1) != nchar(seq2)) {
        stop("seq1 and seq2 must be equal in length")
    }

    # Define property distances
    #dist_mat <- getPropertyMatrix(property)
    if (property == "hydropathy") {
        prop_data <- HYDROPATHY
    } else if (property == "bulkiness") {
        prop_data <- BULKINESS
    } else if (property == "polarity") {
        prop_data <- POLARITY
    }
    
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
    #d <- sum(sapply(1:n, function(x) { dist_mat[seq1[x], seq2[x]] }))
    d <- sum(prop_data[seq2] - prop_data[seq1])
    
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
    if(is.null(hydropathy)) {
        data(aaindex, package="seqinr", envir=environment())
        hydropathy <- aaindex[[151]]$I
    }
    names(hydropathy) <- alakazam::translateStrings(names(hydropathy),
                                                    AA_TRANS)
    # Create character vector from string
    aa <- strsplit(as.character(seq), "")
    # Function to translate a single string
    .gravy <- function(x) {
        # Ignore AAs that are "*" and "X"
        x <- x[x != "X" & x != "*"]
        sum(hydropathy[x])/length(x)
    }
    aa_gravy <- sapply(aa, .gravy)
    return(aa_gravy)
}


#' Calculates aliphatic index of amino acid sequences.
#'
#' @param    seq           vector of amino acid sequences.
#' 
#' @return   aliphatic index for the sequence.
#'
#' @examples 
#' seq <- "CARDRSTPWRRGIASTTVRTSW"
#' aliphatic(seq)
#'
#' @export
aliphatic <- function(seq) {
    # TODO: Why are we normalizing aliphatic index by length?
    n_ala <- countOccurrences(seq, "[A]")
    n_val <- countOccurrences(seq, "[V]")
    n_leu_ile <- countOccurrences(seq, "[LI]")
    a <- 2.9
    b <- 3.9
    aa_aliphatic = (n_ala + a*n_val + b*n_leu_ile) / nchar(seq)
    return(aa_aliphatic)
}


#' Calculates charge of amino acid sequences.
#'
#' @param    seq           vector of amino acid sequences.
#' @param    ph            environmental pH.
#' @param    pk            named vector defining pK values for each amino acid.
#' 
#' @return   charge of the sequence.
#'
#' @examples 
#' seq <- "CARDRSTPWRRGIASTTVRTSW"
#' charge(seq)
#'
#' @export
charge <- function(seq, ph=7, pk=NULL) {
    if(is.null(pk)) {
        data(pK, package="seqinr", envir=environment())
        pk <- pK$EMBOSS
        names(pk) <- rownames(pK)
    }
    names(pk) <- alakazam::translateStrings(names(pk), AA_TRANS)
    carg <- countOccurrences(seq, "R") * (1/(1 + 10^(1 * (ph - pk["R"]))))
    chis <- countOccurrences(seq, "H") * (1/(1 + 10^(1 * (ph - pk["H"]))))
    clys <- countOccurrences(seq, "K") * (1/(1 + 10^(1 * (ph - pk["K"]))))
    casp <- countOccurrences(seq, "D") * (-1/(1 + 10^(-1 * (ph - pk["D"]))))
    cglu <- countOccurrences(seq, "E") * (-1/(1 + 10^(-1 * (ph - pk["E"]))))
    ccys <- countOccurrences(seq, "C") * (-1/(1 + 10^(-1 * (ph - pk["C"]))))
    ctyr <- countOccurrences(seq, "Y") * (-1/(1 + 10^(-1 * (ph - pk["Y"]))))
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
#'                               "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
#' regionProperties(db, sequenceColumn="JUNCTION", nt=TRUE, trim=TRUE, region_name="CDR3")
#' 
#' @export
regionProperties <- function(data, seq="JUNCTION", 
                             nt=FALSE, trim=FALSE, label=NULL) {
    #nt=T
    #trim=T
    #region_name=NULL
    
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