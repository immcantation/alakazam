# Common data structure manipulation, file input/output, and plotting functions for Alakazam
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2014.9.24

#' Read a Change-O tab delimited database file
#' 
#' Reads a tab delimited database file created by a Change-O tool into a 
#' a \code{data.frame}.
#'
#' @param    file       a tab delimited database file output by a Change-O tool.
#' @param    seq_upper  if \code{TRUE} convert sequence fields to upper case.
#'                      if \code{FALSE} do not alter sequence fields.
#' @return   a \code{data.frame} of the database file
#' 
#' @seealso  \code{\link{read.table}}
#' @examples
#' \dontrun{
#'   # Loads a database and converts sequences to uppercase characters
#'   df <- readChangeoDb("changeo_output.tab")
#' 
#'   # Do not alter sequence field case
#'   df <- readChangeoDb("changeo_output.tab", seq_upper=FALSE)
#' }
#' 
#' @export
readChangeoDb <- function(file, seq_upper=TRUE) {
    # Read file
    db_df <- read.delim(file, as.is=TRUE, na.strings=c("", "NA", "None"))
    
    # Handled genotyped data
    # ifelse("V_CALL_GENOTYPED" %in% colnames(data), "V_CALL_GENOTYPED", "V_CALL")
    
    # Convert sequence fields to upper case
    if (seq_upper) {
        seq_columns <- c("SEQUENCE", "SEQUENCE_GAP", "JUNCTION", 
                         "GERMLINE_GAP", "GERMLINE_GAP_D_MASK")
        
        for (x in intersect(seq_columns, colnames(db_df))) {
            db_df[, x] <- toupper(db_df[, x]) 
        }
    }
    
    return(db_df)
}


#' Create a temporary folder
#'
#' \code{makeTempDir} creates a randomly named temporary folder in the 
#' system temp location.
#' 
#' @param    prefix  a prefix name for the folder.
#' @return   The path to the temporary folder.
#' 
#' @seealso  This is just a wrapper for \code{\link{tempfile}} and 
#'           \code{\link{dir.create}}.
#' @examples
#' makeTempDir("Clone50")
#' 
#' @export
makeTempDir <- function(prefix) {
    temp_path <- tempfile(paste0(prefix, "-temp-"))
    dir.create(temp_path)
    
    return(temp_path)
}


#' Translate a vector of strings
#' 
#' \code{translateStrings} modifies a character vector by substituting one or more 
#' strings with a replacement string.
#'
#' Replacement is accomplished using \code{\link{gsub}}. 
#' 
#' Does not perform partial replacements. Each translation value must match a complete 
#' \code{strings} value or it will not be replaced.  Values that do not have a replacement
#' named in the \code{translation} parameter will not be modified.
#' 
#' @param    strings      a vector of strings to modify.
#' @param    translation  a named character vector or a list of character vectors specifying 
#'                        the strings to replace (values) and their replacements (names).
#' @return   A modified \code{strings} vector.
#' 
#' @seealso  \code{\link{gsub}} for single value replacement in the base package.
#' @examples
#' # Using a vector translation
#' strings <- LETTERS[1:5]
#' translation <- c("POSITION1"="A", "POSITION5"="E")
#' translateStrings(strings, translation)
#' 
#' # Using a list translation
#' strings <- LETTERS[1:5]
#' translation <- list("1-3"=c("A","B","C"), "4-5"=c("D","E"))
#' translateStrings(strings, translation)
#' 
#' @export
translateStrings <- function(strings, translation) {
    for (n in names(translation)) {
        rep_regex <- paste(translation[[n]], collapse='|')
        strings <- gsub(paste0("^(", rep_regex, ")$"), n, strings)
    }
    
    return(strings)
}
