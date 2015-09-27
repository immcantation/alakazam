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
# @param    n  a vector of n.
# @param    k  a vector of k.
#
# @return   The approximation of log(n choose k). For n < 100 \link{lchoose} is used.
lchooseStirling <- function(n, k) {
    if (any(n < k)) {
        stop("n must be >= k")
    }
    
    n_len <- length(n)
    k_len <- length(k)
    nCk <- rep(NA, max(n_len, k_len))
    nCk[n == k] <- 0
    
    # a = index n_small
    # i = index k_small
    # x = index nCk_small
    #
    # b = index n_large
    # j = index k_large
    # y = index nCk_large
    #
    # Check for vector inputs and assign indexing
    if (n_len >= 1 & k_len >= 1 & n_len == k_len) {
        a <- i <- x <- (n < 100 & n != k)
        b <- j <- y <- (n >= 100 & n != k)
    } else if (n_len > 1 & k_len == 1) {
        a <- x <- (n < 100 & n != k)
        b <- y <- (n >= 100 & n != k)
        i <- j <- TRUE
    } else if (n_len == 1 & k_len > 1) {
        a <- (n < 100)
        b <- !a
        i <- j <- (n != k)
        x <- if (n < 100) { i } else { NULL }
        y <- if (n >= 100) { i } else { NULL }
    } else {
        stop("Inputs are wrong. n and k must have the same length or be length one.")
    } 

        
    # Small n
    nCk[x] <-  lchoose(n[a], k[i])
    
    # Large n indices
    nCk[y] <- n[b]*log(n[b]) - k[j]*log(k[j]) - (n[b] - k[j])*log(n[b] - k[j]) + 
              0.5*(log(n[b]) - log(k[j]) - log(n[b] - k[j]) - log(2*pi))
    
#     .nCk <- function(n, k) {
#         n*log(n) - k*log(k) - (n - k)*log(n - k) + 
#         0.5*(log(n) - log(k) - log(n - k) - log(2*pi))
#     }
    
    return(nCk)
}