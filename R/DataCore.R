#' Common data structure manipulation, file input/output, and plotting functions for Alakazam
#' 
#' @author     Jason Anthony Vander Heiden
#' @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @version    0.2.0
#' @date       2014.9.18


#### Imports ####
require(grid)

#### Change-O file input/output functions ####

#' Read a Change-O tab delimited database file
#'
#' @param     file         a tab delimited file output by Change-O
#' @param     seq_upper    if TRUE convert sequence fields to upper case
#'                         if FALSE do not alter sequence fields
#' @return    a data.frame of the database file
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


#' Create temporary folder
#'
#' @param     prefix    a prefix name for the folder
#' @return    the path to the temporary folder
makeTempDir <- function(prefix) {
    temp_path <- tempfile(paste0(prefix, "-temp-"))
    dir.create(temp_path)
    
    return(temp_path)
}


#### Data structure manipulation functions ####

#' Transalate a set of strings into single replacement strings
#'
#' @param     strings        a vector of strings to modify
#' @param     translation    a list of  strings vectors specifying the strings 
#'                           to replace (values) and their replacement (names)
#' @return    NULL
translateStrings <- function(strings, translation) {
    for (n in names(translation)) {
        strings <- gsub(paste(translation[[n]], collapse='|'), n, strings)
    }
    
    return(strings)
}


#### Plotting functions ####

#' Plot multiple ggplot objects
#' 
#' @param     ...     ggplot objects to plot
#' @param     ncol    number of columns in the plot 
#' @return    NULL
#' 
#' @references  http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
multiplot <- function(..., ncol=1) {
    p <- list(...)
    n <- length(p)
    layout <- matrix(seq(1, ncol*ceiling(n/ncol)), ncol=ncol, nrow=ceiling(n/ncol))
    
    # Plot
    if (n == 1) {
        plot(p[[1]])
    } else {
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:n) {
            idx <- as.data.frame(which(layout == i, arr.ind=T))
            plot(p[[i]], vp=viewport(layout.pos.row = idx$row, layout.pos.col=idx$col))
        }
    }
}
