# Common input/output and data structure manipulation functions for Alakazam
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2015.03.05


#### File I/O functions ####

#' Read a Change-O tab delimited database file
#' 
#' \code{readChangeoDb} reads a tab delimited database file created by a Change-O tool 
#' into a data.frame.
#'
#' @param    file       tab delimited database file output by a Change-O tool.
#' @param    select     columns to select from database file.
#' @param    drop       columns to drop from database file.
#' @param    seq_upper  if \code{TRUE} convert sequence fields to upper case;
#'                      if \code{FALSE} do not alter sequence fields.
#' 
#' @return   A data.frame of the database file.
#' 
#' @seealso  See \code{\link{writeChangeoDb}} for file output. 
#'           Uses \code{\link{read.table}}.
#' @family   file input and output functions
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' 
#' # Subset columns and convert sequence fields to upper case
#' df <- readChangeoDb(file, select=c("SEQUENCE_ID", "SEQUENCE_GAP"), seq_upper=TRUE)
#' 
#' # Drop columns and do not alter sequence field case
#' df <- readChangeoDb(file, drop=c("D_CALL", "DUPCOUNT"), seq_upper=FALSE)
#' 
#' @export
readChangeoDb <- function(file, select=NULL, drop=NULL, seq_upper=TRUE) {
    # Define column data types
    seq_columns <- c("SEQUENCE", "SEQUENCE_GAP", "JUNCTION", 
                     "GERMLINE_GAP", "GERMLINE_GAP_D_MASK")
    text_columns <- c("SEQUENCE_ID", "CLONE", "SAMPLE")
    
    # Read file
    db_df <- read.delim(file, as.is=TRUE, na.strings=c("", "NA", "None"))
    
    # Select columns
    select_columns <- colnames(db_df)
    if(!is.null(select)) { select_columns <- intersect(select_columns, select) }
    if(!is.null(drop)) { select_columns <- setdiff(select_columns, drop) }
    db_df <- subset(db_df, select=select_columns)
    
    # Convert sequence fields to upper case
    if (seq_upper) {
        for (x in intersect(seq_columns, select_columns)) {
            db_df[, x] <- toupper(db_df[, x]) 
        }
    }
    
    # Convert text fields to character
    for (x in intersect(text_columns, select_columns)) {
        db_df[, x] <- as.character(db_df[, x])
    }
    
    return(db_df)
}


#' Write a Change-O tab delimited database file
#' 
#' \code{writeChangeoDb} is a simple wrapper around \code{\link{write.table}} with defaults 
#' appropriate for writing a Change-O tab delimited database file from a data.frame.
#'
#' @param    data  data.frame of Change-O data.
#' @param    file  output file name.
#' 
#' @return   NULL
#' 
#' @seealso  See \code{\link{readChangeoDb}} for file input. 
#'           Uses \code{\link{write.table}}.
#' @family   file input and output functions
#' 
#' @examples
#' \dontrun{
#'   # Write a database
#'   writeChangeoDb(data, "changeo_output.tab")
#' }
#' 
#' @export
writeChangeoDb <- function(data, file) {
    write.table(data, file=file, quote=FALSE, sep="\t", row.names=FALSE)
}


#' Create a temporary folder
#'
#' \code{makeTempDir} creates a randomly named temporary folder in the 
#' system temp location.
#' 
#' @param    prefix  prefix name for the folder.
#' 
#' @return   The path to the temporary folder.
#' 
#' @seealso  This is just a wrapper for \code{\link{tempfile}} and 
#'           \code{\link{dir.create}}.
#' @family   file input and output functions
#' 
#' @examples
#' makeTempDir("Clone50")
#' 
#' @export
makeTempDir <- function(prefix) {
    temp_path <- tempfile(paste0(prefix, "-temp-"))
    dir.create(temp_path)
    
    return(temp_path)
}


#### Data manipulation functions ####

#' Translate a vector of strings
#' 
#' \code{translateStrings} modifies a character vector by substituting one or more 
#' strings with a replacement string.
#'
#' @param    strings      vector of character strings to modify.
#' @param    translation  named character vector or a list of character vectors specifying 
#'                        the strings to replace (values) and their replacements (names).
#' 
#' @return   A modified \code{strings} vector.
#' 
#' @details
#' Does not perform partial replacements. Each translation value must match a complete 
#' \code{strings} value or it will not be replaced.  Values that do not have a replacement
#' named in the \code{translation} parameter will not be modified.
#' 
#' Replacement is accomplished using \code{\link{gsub}}.
#' 
#' @seealso  See \code{\link{gsub}} for single value replacement in the base package.
#' 
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
    # TODO:  use match instead for complete matching?  Currently regex characters in values will mess up the matching.
    for (n in names(translation)) {
        rep_regex <- paste(translation[[n]], collapse='|')
        strings <- gsub(paste0("^(", rep_regex, ")$"), n, strings)
    }
    
    return(strings)
}
