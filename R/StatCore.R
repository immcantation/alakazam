# Common statistics functions for Alakazam

#' Weighted meta-analysis of p-values via Stouffer's method
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

# Stirling's approximation of the binomial coefficient
# 
# Calculates Stirling's approximation of the binomial coefficient for large numbers.
#
# @param    n  n.
# @param    k  k.
#
# @return   The approximation of log(n choose k).
lchooseStirling <- function(n, k) {
    x <- n*log(n) - k*log(k) - (n - k)*log(n - k) + 
         0.5*(log(n) - log(k) - log(n - k) - log(2*pi))
    
    return(x)
}