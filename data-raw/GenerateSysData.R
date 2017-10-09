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

# Changeo data type specification
CHANGEO <- list(SEQUENCE_ID="c",
                SEQUENCE_INPUT="c",
                FUNCTIONAL="l",
                IN_FRAME="l",
                STOP="l",
                MUTATED_INVARIANT="c",
                INDELS="l",
                V_CALL="c",
                D_CALL="c",
                J_CALL="c",
                SEQUENCE_VDJ="c",
                SEQUENCE_IMGT="c",
                V_SEQ_START="i",
                V_SEQ_LENGTH="i",
                V_GERM_START_VDJ="i",
                V_GERM_LENGTH_VDJ="i",
                V_GERM_START_IMGT="i",
                V_GERM_LENGTH_IMGT="i",
                NP1_LENGTH="i",
                D_SEQ_START="i",
                D_SEQ_LENGTH="i",
                D_GERM_START="i",
                D_GERM_LENGTH="i",
                NP2_LENGTH="i",
                J_SEQ_START="i",
                J_SEQ_LENGTH="i",
                J_GERM_START="i",
                J_GERM_LENGTH="i",
                JUNCTION_START="i",
                JUNCTION_LENGTH="i",
                JUNCTION="c",
                JUNCTION_AA="c",
                FWR1_IMGT="c",
                FWR2_IMGT="c",
                FWR3_IMGT="c",
                FWR4_IMGT="c",
                CDR1_IMGT="c",
                CDR2_IMGT="c",
                CDR3_IMGT="c",
                V_SCORE="d",
                V_IDENTITY="d",
                V_EVALUE="d",
                V_BTOP="c",
                D_SCORE="d",
                D_IDENTITY="d",
                D_EVALUE="d",
                D_BTOP="c",
                J_SCORE="d",
                J_IDENTITY="d",
                J_EVALUE="d",
                J_BTOP="c",
                HMM_SCORE="d",
                N1_LENGTH="i",
                N2_LENGTH="i",
                P3V_LENGTH="i",
                P5D_LENGTH="i",
                P3D_LENGTH="i",
                P5J_LENGTH="i",
                D_FRAME="i",
                CDR3_IGBLAST_NT="c",
                CDR3_IGBLAST_AA="c",
                GERMLINE_VDJ="c",
                GERMLINE_VDJ_V_REGION="c",
                GERMLINE_VDJ_D_MASK="c",
                GERMLINE_IMGT="c",
                GERMLINE_IMGT_V_REGION="c",
                GERMLINE_IMGT_D_MASK="c",
                GERMLINE_V_CALL="c",
                GERMLINE_D_CALL="c",
                GERMLINE_J_CALL="c",
                GERMLINE_REGIONS="c",
                V_CALL_GENOTYPED="c",
                CLONE="c",
                PRIMER="c",
                PRCONS="c",
                CONSCOUNT="i",
                DUPCOUNT="i")


#### Save to R/sysdata.rda ####

devtools::use_data(HYDROPATHY_KYTJ82,
                   BULKINESS_ZIMJ68,
                   POLARITY_GRAR74,
                   PK_EMBOSS,
                   CHANGEO,
                   internal=TRUE, overwrite=TRUE)