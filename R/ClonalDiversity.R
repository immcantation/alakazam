#' Clonal diversity analysis
#' 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2014.9.1


#### Class definitions ####

#' S4 class defining diversity curve over multiple diversity orders
setClass("DiversityCurve", contains="data.frame",
         slots=c(groups="character", n="numeric", nboot="numeric", ci="numeric"))

#' S4 class defining diversity signficance test at fixed diversity order
setClass("DiversityTest", contains="data.frame",
         slots=c(groups="character", stats="data.frame", n="numeric", 
                 nboot="numeric", q="numeric"))

#### Calculation functions ####

#' Calculate the diversity index
#'
#' @param         p    a numerical vector of species counts or proportions
#' @param         q    a numerical vector of diversity orders
#' @return        a vector of diversity numbers (D) for each q
#' 
#' @references    Hill, M. Diversity and evenness: a unifying notation and its consequences. 
#'                    Ecology 54, 427–432 (1973).
calcDiversity <- function(p, q) {
    # Convert p to proportional abundance
    p <- p / sum(p)
    # Add jitter to q=1
    q[q == 1] <- 0.9999
    # Calculate D for each q
    D <- sapply(q, function(x) sum(p^x)^(1 / (1 - x)))
    
    return(D)
}


#' Calculates the diversity index curve for clones
#'
#' @param     data      a data.frame containing clonal assignments
#' @param     group     the data column containing group identifiers
#' @param     clone     the data column containing clone identifiers
#' @param     min_q     the minimum value of q
#' @param     max_q     the maximum value of q
#' @param     step_q    the value to increment q
#' @param     min_n     the minimum number of observations to sample
#'                      a group with less observations than the minimum is excluded
#' @param     max_n     the maximum number of observations to sample
#' @param     ci        the confidence interval to calculate
#' @param     nboot     the number of bootstrap iterations to perform
#' @return    a data.frame summarizing the diversity scores over all q
bootstrapDiversity <- function(data, group, clone="CLONE", min_q=0, max_q=32, step_q=0.05, 
                               min_n=10, max_n=NULL, ci=0.95, nboot=2000) {
    # >>> Add MAD
    
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
    cat("-> CALCULATING DIVERSITY\n")
    pb <- txtProgressBar(min=1, max=length(group_keep), initial=1, width=40, style=3)
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


#' Performs pairwise tests of the diversity index at a given diversity order
#'
#' @param     data      a data.frame containing clonal assignments
#' @param     q         the diversity order to test
#' @param     group     the data column containing group identifiers
#' @param     clone     the data column containing clone identifiers
#' @param     min_n     the minimum number of observations to sample
#'                      a group with less observations than the minimum is excluded
#' @param     max_n     the maximum number of observations to sample
#' @param     nboot     the number of bootstrap iterations to perform
#' @return    a data.frame summarizing the diversity scores over all q
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

#' Define universal plot settings
#'
#' @return    a ggplot2 theme object
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


#' Plots the results of bootstrapDiversity
#'
#' @param     data            a DiversityCurve object returned by boostrapDiversity
#' @param     colors          a named character vector whose names are values in the group column
#'                            and whose values are colors to assign to those group values
#' @param     main_title      the plot title
#' @param     legend_title    the legend title          
#' @param     log_q           if TRUE then plot q in a log scale
#'                            if FALSE plot in a linear scale
#' @param     log_d           if TRUE then plot diversity in a log scale
#'                            if FALSE plot in a linear scale
#' @return    NULL                
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
}

