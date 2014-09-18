#' Common DNA, amino acid, and gene annotation operations for Alakazam
#' 
#' @author     Jason Anthony Vander Heiden
#' @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @version    0.2.0
#' @date       2014.9.17

#### Imports ####
#require(seqinr)
require(plyr)
require(stringr)

#### Constants ####

# Default colors
NUC_COLORS <- c("A"="#64F73F", 
                "C"="#FFB340", 
                "G"="#EB413C", 
                "T"="#3C88EE",
                "U"="#3C88EE")
IG_COLORS <- c("IgA"="#377EB8", 
               "IgD"="#FF7F00", 
               "IgE"="#E41A1C", 
               "IgG"="#4DAF4A", 
               "IgM"="#984EA3",
               "IgK"="#E5C494",
               "IgL"="#FFD92F")
TR_COLORS <- c("TRA"="#CBD5E8", 
               "TRB"="#F4CAE4", 
               "TRD"="#FDCDAC", 
               "TRG"="#E6F5C9")


#### DNA functions ####

#### Amino acid functions ####

#### Gene annotation functions ####

#' Get segment family calls
#'
#' @param     segment_call    a character vector containing segment calls delimited by commas
#' @param     first           if TRUE return only the first call
#'                            if FALSE return all calls delimited by commas
#' @return    a string containg the family calls
getFamily <- function(segment_call, first=TRUE) {
    family_regex <- '(IG[HLK][VDJ]\\d+)'
    
    if (first) {
        r <- str_extract(segment_call, perl(family_regex))
    } else {
        r <- str_extract_all(segment_call, perl(family_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}


#' Get segment gene calls
#'
#' @param     segment_call    a character vector containing segment calls delimited by commas
#' @param     first           if TRUE return only the first call
#'                            if FALSE return all calls delimited by commas
#' @return    a string containg the gene calls
getGene <- function(segment_call, first=TRUE) {
    gene_regex <- '(IG[HLK][VDJ]\\d+[-/\\w]*)'
    
    if (first) {
        r <- str_extract(segment_call, perl(gene_regex))
    } else {
        r <- str_extract_all(segment_call, perl(gene_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}


#' Get segment allele calls
#'
#' @param     segment_call    a string containing one or more segment calls delimited by commas
#' @param     first           if TRUE return only the first call
#'                            if FALSE return all calls delimited by commas
#' @return    a string containg the allele calls
getAllele <- function(segment_call, first=TRUE) {
    allele_regex <- '(IG[HLK][VDJ]\\d+[-/\\w]*[-\\*][\\.\\w]+)'
    
    if (first) {
        r <- str_extract(segment_call, perl(allele_regex))
    } else {
        r <- str_extract_all(segment_call, perl(allele_regex))
        r <- sapply(r, function(x) paste(sort(x), collapse=','))
    }
    
    return(r)
}
