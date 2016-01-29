# Generate sysdata

#### Default amino acid property set ####

# Wrapper function to pull amino acid property vectors from seqinr
#
# @param    property  string defining the property. One of
#                     \code(c("hydropathy", "bulkiness", "polarity", "pK"))
#
# @return   A vector of scores for the given property with single
#           character amino acid labels.
#
# @examples
# getPropertyData("hydro")
# getPropertyData("bulk")
# getPropertyData("polar")
# getPropertyData("pK")
getPropertyData <- function(property){
    property <- match.arg(property, c("hydropathy", "bulkiness", "polarity" ,"pK"))
    
    # Setup new environment to avoid R CMD check NOTE
    e1 <- new.env(parent=environment())
    data(aaindex, package="seqinr", envir=e1)
    data(pK, package="seqinr", envir=e1)
    
    if (property == "hydropathy") {
        # Kyte & Doolittle, 1982.
        scores <- with(e1, aaindex[["KYTJ820101"]]$I)
        names(scores) <- translateStrings(names(scores), ABBREV_AA)
    } else if (property == "bulkiness") {
        # Zimmerman et al, 1968.
        scores <- with(e1, aaindex[["ZIMJ680102"]]$I)
        names(scores) <- translateStrings(names(scores), ABBREV_AA)
    } else if (property == "polarity") {
        # Grantham, 1974
        scores <- with(e1, aaindex[["GRAR740102"]]$I)
        names(scores) <- translateStrings(names(scores), ABBREV_AA)
    } else if (property == "pK") {
        # EMBOSS
        scores <- with(e1, setNames(pK[["EMBOSS"]], rownames(pK)))
    } 
    
    return(scores)
}

# Assign to internal variables
HYDROPATHY_KYTJ82 <- getPropertyData("hydropathy")
BULKINESS_ZIMJ68 <- getPropertyData("bulkiness")
POLARITY_GRAR74 <- getPropertyData("polarity")
PK_EMBOSS <- getPropertyData("pK")

#### Save to R/sysdata.rda ####

devtools::use_data(HYDROPATHY_KYTJ82,
                   BULKINESS_ZIMJ68,
                   POLARITY_GRAR74,
                   PK_EMBOSS,
                   internal=TRUE, overwrite=TRUE)