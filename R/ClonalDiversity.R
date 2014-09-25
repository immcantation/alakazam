# Clonal diversity analysis
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2014.9.25


#### Classes ####

#' S4 class defining diversity curve 
#'
#' \code{DiversityCurve} defines diversity (D) scores over multiple diversity orders (q).
#' 
#' @slot  .Data   data.frame with columns (q, D, lower, upper) defining the median (D),
#'                upper confidence bound (upper) lower confidence bound (lower) for each
#'                value of q (q).
#' @slot  groups  character vector of groups retained in diversity calculation.
#' @slot  n       numeric value indication the number of sampled sequences from each group.
#' @slot  nboot   number of bootstrap realizations.
#' @slot  ci      confidence interval (between 0 and 1).
#' 
#' @name DiversityCurve
#' @export
setClass("DiversityCurve", contains="data.frame",
         slots=c(groups="character", n="numeric", nboot="numeric", ci="numeric"))


#' S4 class defining diversity significance
#'
#' \code{DiversityTest} defines the signifance of diversity (D) differences at a fixed
#' fixed diversity order (q).
#' 
#' @slot  .Data   data.frame with columns for group pairs and associated pvalues, and 
#'                bootstrap delta distribution summary statistics (delta_mean, delta_sd).
#' @slot  groups  character vector of groups retained in diversity calculation.
#' @slot  stats   data.frame of diversity (D) summary statistics by group. Includes columns
#'                for group, mean, median, sd, mad.
#' @slot  n       numeric value indication the number of sampled sequences from each group.
#' @slot  nboot   number of bootstrap realizations.
#' @slot  q       diversity order tested (q).
#' 
#' @name DiversityTest
#' @export
setClass("DiversityTest", contains="data.frame",
         slots=c(groups="character", stats="data.frame", n="numeric", 
                 nboot="numeric", q="numeric"))

#### Calculation functions ####

#' Calculate the diversity index
#' 
#' \code{calcDiversity} calculates the diversity index for a vector of diversity orders.
#'
#' @param    p  numeric vector of species counts or proportions.
#' @param    q  numeric vector of diversity orders.
#' @return   A vector of diversity numbers (D) for each q.
#' 
#' @seealso  Used by \code{\link{bootstrapDiversity}} and \code{\link{testDiversity}}.
#' @references    
#' Hill, M. Diversity and evenness: a unifying notation and its consequences. 
#'   Ecology 54, 427–432 (1973).
#' @examples
#' # May define p as clonal member counts
#' p <- c(1, 1, 3, 10)
#' q <- c(0, 1, 2)
#' calcDiversity(p, q)
#'
#' # Or proportional abundance
#' p <- c(1/15, 1/15, 1/5, 2/3)
#' calcDiversity(p, q)
#' 
#' @export
calcDiversity <- function(p, q) {
    # Convert p to proportional abundance
    p <- p / sum(p)
    # Add jitter to q=1
    q[q == 1] <- 0.9999
    # Calculate D for each q
    D <- sapply(q, function(x) sum(p^x)^(1 / (1 - x)))
    
    return(D)
}


#' Generate a clonal diversity index curve
#'
#' \code{bootstrapDiversity} divides a set of clones by a group annotation,
#' uniformly subsamples the sequences from each group, and calculates diversity
#' scores (D) over an interval of diversity orders (q).
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    min_q     minimum value of q.
#' @param    max_q     maximum value of q.
#' @param    step_q    value to increment q.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} the maximum
#'                     if automatically determined from the size of the largest group.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap iterations to perform.
#' @return   A \code{DiversityCurve} summarizing the diversity scores over all q.
#' 
#' @seealso  See \code{\link{calcDiversity}} for the basic calculation and 
#'           \code{\link{DiversityCurve}} for the return object. 
#'           See \code{\link{testDiversity}} for significance testing.
#' @examples
#' # Load example data
#' file <- system.file("extdata", "IB_T_genotyped_clone-pass_germ-pass_300.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # All groups pass default minimum sampling threshold of 10 sequences
#' div <- bootstrapDiversity(df, "BARCODE", step_q=1, max_q=10, nboot=100)
#' div[div$q == 0, ]
#' slot(div, "groups")
#' slot(div, "n")
#' 
#' # Increasing threshold results in exclusion of small groups and a warning message
#' div <- bootstrapDiversity(df, "BARCODE", min_n=40, step_q=1, max_q=10, nboot=100)
#' div[div$q == 0, ]
#' slot(div, "groups")
#' slot(div, "n")
#'
#' @export
bootstrapDiversity <- function(data, group, clone="CLONE", min_q=0, max_q=32, step_q=0.05, 
                               min_n=10, max_n=NULL, ci=0.95, nboot=2000) {
    # Verify function arguments
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    if (!(group %in% names(data))) {
        stop(paste("The column", group, "does not exist in the input data.frame"))
    }
    if (!(clone %in% names(data))) {
        stop(paste("The column", clone, "does not exist in the input data.frame"))
    }
    
    # Count observations per group and set sampling criteria
    group_sum <- ddply(data, c(group), here(summarize), count=length(eval(parse(text=group))))
    group_all <- as.character(group_sum[, group])
    group_keep <- as.character(group_sum[group_sum$count >= min_n, group])
    n <- min(group_sum$count[group_sum$count >= min_n], max_n)
    q <- seq(min_q, max_q, step_q)
    ci_probs <- c((1 - ci)/2, 0.5, ci + (1 - ci)/2)

    # Warn if groups removed
    if (length(group_keep) < length(group_all)) {
        warning("Not all groups passed min_n=", min_n, " threshold. Excluded: ", 
                paste(setdiff(group_all, group_keep), collapse=", "))
    }
    
    # Generate diversity index and confidence intervals via resampling
    # >>> Add MAD?
    cat("-> CALCULATING DIVERSITY\n")
    pb <- txtProgressBar(min=0, max=length(group_keep), initial=0, width=40, style=3)
    div_list <- list()
    i <- 0
    for (g in group_keep) {
        i <- i + 1
        r <- which(data[[group]] == g)
        sample_mat <- replicate(nboot, data[[clone]][sample(r, n, replace=T)])
        boot_mat <- apply(sample_mat, 2, function(x) calcDiversity(table(x), q))
        boot_ci <- t(apply(boot_mat, 1, quantile, probs=ci_probs))
        div_list[[g]] <- matrix(c(q, boot_ci[, c(2, 1, 3)]),
                                nrow=length(q), ncol=4, 
                                dimnames=list(NULL, c("q", "D", "lower", "upper")))
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    
    # Generate return object
    div <- new("DiversityCurve", ldply(div_list, .id="group"), 
               groups=group_keep, n=n, nboot=nboot, ci=ci)

    return(div)
}


#' Pairwise test of the diversity index
#' 
#' \code{testDiversity} performs pairwise significance tests of the diversity index (D)
#' at a given diversity order (q) for a set of annotation groups using a bootstrap delta
#' distribution.
#'
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    q         diversity order to test.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} the maximum
#'                     if automatically determined from the size of the largest group.
#' @param    nboot     number of bootstrap realizations to perform.
#' @return   A \code{DiversityTest} object containing p-values and summary statistics.
#' 
#' @seealso  See \code{\link{calcDiversity}} for the basic calculation and 
#'           \code{\link{DiversityTest}} for the return object. 
#'           See \code{\link{bootstrapDiversity}} for curve generation.
#' @examples          
#' # Load example data
#' file <- system.file("extdata", "IB_T_genotyped_clone-pass_germ-pass_300.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Groups under the size threshold are excluded and a warning message is issued.
#' div <- testDiversity(df, "BARCODE", q=0, min_n=30, nboot=100)
#' div
#' slot(div, "groups")
#' slot(div, "n")
#' slot(div, "stats")
#' 
#' @export
testDiversity <- function(data, q, group, clone="CLONE", min_n=10, max_n=NULL, nboot=2000) {
    
    #>>> Needs two-tailed variant.  P-values can >1 if pvalue=pvalue*2
    
    # Verify function arguments
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    if (!(group %in% names(data))) {
        stop(paste("The column", group, "does not exist in the input data.frame"))
    }
    if (!(clone %in% names(data))) {
        stop(paste("The column", clone, "does not exist in the input data.frame"))
    }
    
    # Count observations per group and set sampling criteria
    group_sum <- ddply(data, c(group), here(summarize), count=length(eval(parse(text=group))))
    group_all <- as.character(group_sum[, group])
    group_keep <- as.character(group_sum[group_sum$count >= min_n, group])
    n <- min(group_sum$count[group_sum$count >= min_n], max_n)
    ngroup <- length(group_keep)
    
    # Warn if groups removed
    if (ngroup < length(group_all)) {
        warning("Not all groups passed min_n=", min_n, " threshold. Excluded: ", 
                paste(setdiff(group_all, group_keep), collapse=", "))
    }
    
    # Generate diversity index and confidence intervals via resampling
    cat("-> CALCULATING DIVERSITY\n")
    pb <- txtProgressBar(min=1, max=ngroup, initial=1, width=40, style=3)
    div_mat <- matrix(NA, nboot, ngroup, dimnames=list(NULL, group_keep))
    for (i in 1:ngroup) {
        g <- group_keep[i]
        r <- which(data[[group]] == g)
        sample_mat <- replicate(nboot, data[[clone]][sample(r, n, replace=T)])
        div_mat[, i] <- apply(sample_mat, 2, function(x) calcDiversity(table(x), q))
        setTxtProgressBar(pb, i)
    }
    cat("\n")
        
    # Compute ECDF of bootstrap distribution shift from bootstrap deltas
    #>>> Change to median and mad?
    group_pairs <- combn(group_keep, 2, simplify=F)
    npairs <- length(group_pairs)
    delta_mat <- matrix(NA, nboot, npairs)
    pvalue_mat <- matrix(NA, npairs, 3, 
                         dimnames=list(NULL, c("pvalue", "delta_mean", "delta_sd")))
    test_names <- character(length=npairs)
    for (i in 1:npairs) {
        g1 <- group_pairs[[i]][1]
        g2 <- group_pairs[[i]][2]
        if (median(div_mat[, g1]) > median(div_mat[, g2])) {
            g_delta <- div_mat[, g1] - div_mat[, g2]
            g_name <- paste(g1, g2, sep=' > ')
        } else {
            g_delta <- div_mat[, g2] - div_mat[, g1]
            g_name <- paste(g2, g1, sep=' > ')
        }  
        test_names[i] <- g_name
        
        # Determine p-value
        g_cdf <- ecdf(g_delta)
        #pvalue_mat[i, ] <- c(g_cdf(0) * 2, mean(g_delta), sd(g_delta))
        pvalue_mat[i, ] <- c(g_cdf(0), mean(g_delta), sd(g_delta))
    }
    
    # >>> Convert mean/sd to slots of DiversityTest?
    # >>> Plot function for sd/mean (dnorm)?
    test_df <- cbind(data.frame(test=test_names), as.data.frame(pvalue_mat))
    stats_df <- data.frame(group=group_keep, 
                           median=apply(div_mat, 2, median),
                           mad=apply(div_mat, 2, mad),
                           mean=apply(div_mat, 2, mean),
                           sd=apply(div_mat, 2, sd))
    
    # Generate return object
    div <- new("DiversityTest", test_df, groups=group_keep, stats=stats_df, 
               n=n, nboot=nboot, q=q)
    
    return(div)
}


#### Plotting functions ####

# Define universal plot settings
#
# @return    a ggplot2 theme object
getBaseTheme <- function() {
    # Define universal plot settings
    base_theme <- theme_bw() + 
        theme(plot.title=element_text(size=16)) +
        theme(text=element_text(size=14)) +
        theme(strip.background=element_rect(fill='white')) + 
        theme(strip.text=element_text(size=16, face='bold')) +
        theme(axis.title=element_text(size=16, vjust=0.25)) +
        theme(axis.text.x=element_text(size=14, vjust=0.5, hjust=0.5)) +
        theme(axis.text.y=element_text(size=14))
    
    return(base_theme)
}


#' Plot the results of bootstrapDiversity
#' 
#' \code{plotDiversityCurve} plots a \code{DiversityCurve} object.
#'
#' @param    data            \code{DiversityCurve} object returned by \code{boostrapDiversity}
#' @param    colors          named character vector whose names are values in the group column
#'                           and whose values are colors to assign to those group values.
#' @param    main_title      string specifying the plot title.
#' @param    legend_title    string specifying the legend title.
#' @param    log_q           if TRUE then plot q in a log scale;
#'                           if FALSE plot in a linear scale.
#' @param    log_d           if TRUE then plot diversity in a log scale;
#'                           if FALSE plot in a linear scale.
#' @return   a \code{ggplot} object
#' 
#' @seealso  Plotting is performed with \code{\link{ggplot}}.
#' @examples
#' # Load example data
#' file <- system.file("extdata", "IB_T_genotyped_clone-pass_germ-pass_300.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # All groups pass default minimum sampling threshold of 10 sequences
#' div <- bootstrapDiversity(df, "BARCODE", min_n=30, step_q=1, max_q=10, nboot=100)
#' plotDiversityCurve(div, main_title=paste(slot(div, "n"), "sequences sampled"),
#'                    legend_title="Barcode")
#' 
#' @export
plotDiversityCurve <- function(data, colors=NULL, main_title="Diversity", 
                               legend_title=NULL, log_q=TRUE, log_d=TRUE) {
    # Define plot elements
    p1 <- ggplot(data, aes(x=q, y=D, group=group)) + 
        ggtitle(main_title) + 
        getBaseTheme() + 
        xlab('q') +
        ylab(expression(''^q * D)) +
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=0.25) +
        geom_line(aes(color=group)) +
        guides(color=guide_legend(title=legend_title), 
               fill=guide_legend(title=legend_title))  
    if (!is.null(colors)) {
        p1 <- p1 + scale_color_manual(name=legend_title, values=colors) +
            scale_fill_manual(name=legend_title, values=colors)
    }        
    if (log_q) {
        p1 <- p1 + scale_x_continuous(trans=log2_trans(),
                                      breaks=trans_breaks('log2', function(x) 2^x),
                                      labels=trans_format('log2', math_format(2^.x)))
    }
    if (log_d) {
        p1 <- p1 + scale_y_continuous(trans=log2_trans(),
                                      breaks=trans_breaks('log2', function(x) 2^x),
                                      labels=trans_format('log2', math_format(2^.x)))
    }
    
    # Plot
    plot(p1)
    
    return(p1)
}

