# Common statistics functions for Alakazam
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2015.03.05


#' Weighted meta analysis of p-values via Stouffer's method
#'
#' \code{stoufferMeta} combines multiple weighted p-values into a meta-analysis p-value
#' using Stouffer's Z-score method.
#' 
#' @param    p   numeric vector of p-values.
#' @param    w   numeric vector of weights.
#' 
#' @return   A named numeric vector with the combined Z-score and p-value in the form
#'           \code{c(Z, pvalue)}.
#' 
#' @examples
#' # Define p-value and weight vectors
#' p <- c(0.1, 0.05, 0.3)
#' w <- c(5, 10, 1)
#'
#' # Unweighted
#' stoufferMeta(p)
#' 
#' # Weighted
#' stoufferMeta(p, w)
#' 
#' @export
stoufferMeta <- function(p, w=NULL) {
    if (is.null(w)) {
        w <- rep(1, length(p))/length(p)
    } else {
        if (length(w) != length(p)) { stop("Length of p and w must equal.") }
        w <- w/sum(w)
    }
    x <- qnorm(1 - p)
    Z  <- sum(w*x) / sqrt(sum(w^2))
    pvalue <- 1 - pnorm(Z)
    return(c(Z=Z, pvalue=pvalue))
}