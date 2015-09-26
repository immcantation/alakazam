# Clonal diversity analysis

#### Classes ####

#' S4 class defining diversity curve 
#'
#' \code{DiversityCurve} defines diversity (\eqn{D}) scores over multiple diversity 
#' orders (\eqn{q}).
#' 
#' @slot  data      data.frame defining the diversity curve with the following columns:
#'                  \itemize{
#'                    \item  \code{q}:      diversity order.
#'                    \item  \code{D}:      median diversity index over all bootstrap 
#'                                          realizations.
#'                    \item  \code{lower}:  lower confidence inverval bound.
#'                    \item  \code{upper}:  upper confidence interval bound.
#'                  }
#' @slot  groups    character vector of groups retained in the diversity calculation.
#' @slot  n         numeric vector indication the number of sequences sampled from each group.
#' @slot  coverage  numeric vector indication the sample coverage for each sampled group.
#' @slot  nboot     number of bootstrap realizations performed.
#' @slot  ci        confidence interval defining the upper and lower bounds 
#'                  (a value between 0 and 1).
#' 
#' @name         DiversityCurve-class
#' @rdname       DiversityCurve-class
#' @aliases      DiversityCurve
#' @exportClass  DiversityCurve
setClass("DiversityCurve", 
         slots=c(data="data.frame", 
                 groups="character", 
                 n="numeric", 
                 coverage="numeric",
                 nboot="numeric", 
                 ci="numeric"))

#' S4 class defining diversity significance
#'
#' \code{DiversityTest} defines the signifance of diversity (\eqn{D}) differences at a 
#' fixed diversity order (\eqn{q}).
#' 
#' @slot  tests    data.frame describing the significance test results with columns:
#'                 \itemize{
#'                   \item  \code{test}:          string listing the two groups tested.
#'                   \item  \code{pvalue}:        p-value for the test.
#'                   \item  \code{delta_median}:  median of the \eqn{D} bootstrap delta 
#'                                                distribution for the test.
#'                   \item  \code{delta_mad}:     median absolute deviation of the \eqn{D} 
#'                                                bootstrap delta distribution for the test.
#'                   \item  \code{delta_mean}:    mean of the \eqn{D} bootstrap delta 
#'                                                distribution for the test.
#'                   \item  \code{delta_sd}:      standard deviation of the \eqn{D} 
#'                                                bootstrap delta distribution for the test.
#'                 }
#' @slot  summary  data.frame containing summary statistics for the diversity index 
#'                 bootstrap distributions, at the given value of \eqn{q}, with columns:
#'                 \itemize{
#'                   \item  \code{group}:   the name of the group.
#'                   \item  \code{median}:  median of the \eqn{D} bootstrap distribution.
#'                   \item  \code{mad}:     median absolute deviation of the \eqn{D} 
#'                                          bootstrap distribution.
#'                   \item  \code{mean}:    mean of the \eqn{D} bootstrap distribution.
#'                   \item  \code{sd}:      standard deviation of the \eqn{D} bootstrap 
#'                                          distribution.
#'                 }
#' @slot  groups   character vector of groups retained in diversity calculation.
#' @slot  q        diversity order tested (\eqn{q}).
#' @slot  n        numeric value indication the number of sequences sampled from each group.
#' @slot  nboot    number of bootstrap realizations.
#' 
#' @name         DiversityTest-class
#' @rdname       DiversityTest-class
#' @aliases      DiversityTest
#' @exportClass  DiversityTest
DiversityTest <- setClass("DiversityTest", 
         slots=c(tests="data.frame",
                 summary="data.frame",
                 groups="character", 
                 q="numeric",
                 n="numeric", 
                 nboot="numeric"))


#### Methods ####

# TODO:  plot method for DiversityCurve pointing to plotDiversityCurve
# TODO:  plot method for DiversityTest
# TODO:  summary method for DiversityTest
# TODO:  summary method for DiversityCurve

#' @param    x  DiversityCurve object
#' 
#' @rdname   DiversityCurve-class
#' @aliases  print,DiversityCurve-method
setMethod("print", "DiversityCurve", function(x) { print(x@data) })

#' @param    x  DiversityTest object
#' 
#' @rdname   DiversityTest-class
#' @aliases  print,DiversityTest-method
setMethod("print", "DiversityTest", function(x) { print(x@tests) })

# @rdname DiversityCurve
# @export
# setMethod("plot", 
#           signature("DiversityCurve", 
#                     colors="character", 
#                     main_title="character", 
#                     legend_title="character", 
#                     log_q="logical", 
#                     log_d="logical",
#                     xlim="numeric", 
#                     ylim="numeric", 
#                     silent="logical"),
#           function(data, colors=NULL, main_title="Diversity", 
#                    legend_title=NULL, log_q=TRUE, log_d=TRUE,
#                    xlim=NULL, ylim=NULL, silent=FALSE) {
#             plotDiversityCurve(data, colors=colors, main_title=main_title, 
#                                legend_title=legend_title, log_q=log_q, log_d=log_d,
#                                xlim=xlim, ylim=ylim, silent=silent) })


#### Calculation functions ####

# Calculate undetected species
# 
# Calculates the lower bound of undetected species counts using the Chao1 estimator.
#
# @param    x  vector of observed abundance counts.
# 
# @return   The count of undetected species.
inferUnseenCount <- function(x) {
    x <- x[x > 0]
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    
    if (f2 > 0) {
        f0 <- ceiling(((n - 1) * f1^2) / (n * 2 * f2))
    } else {
        f0 <- ceiling(((n - 1) * f1 * (f1 - 1)) / (n * 2))
    }
    
    return(f0)
}


# Calculates diversity under rarefaction
# 
# Calculates Hill numbers under rarefaction
#
# @param    x  vector of observed abundance counts.
# @param    q  numeric vector of diversity orders.
# @param    n  the sequence count to rarefy to.
#
# @return   A vector of diversity scores \eqn{D} for each \eqn{q}.
calcRarefiedDiversity <- function(x, q, n) {
    x <- x[x >= 1]
    N <- sum(x)
    if (n > N) {
        stop("n must be <= the total count of observed sequences.")
    }
    q[q == 1] <- 0.9999
    
    # TODO:  Walk through the math carefully and make sure this is correct!!!
    fk <- sapply(1:n, function(k) sum(exp(lchoose(x, k) + lchoose(N - x, n - k) - lchoose(N, n))))
    D <- sapply(q, function(z) sum((1:n / n)^z * fk)^(1 / (1 - z)))
    
    return(D)
}

# Calculate sample coverage
# 
# Calculates sample coverage of varying orders.
#
# @param    x  numeric vector of abundance as counts
# @param    r  the coverage order to calculate
# 
# @return   The sample coverage of the given order \code{r}.
calcCoverage <- function(x, r=1) {
    n <- sum(x)
    fr <- sum(x == r)
    fs <- sum(x == r + 1)
    
    if (fr == 0) {
        stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", r, ".")
    }
    if (fs == 0) {
        stop("Cannot calculate coverage of order ", r, ". No abundance data with count=", r + 1, ".")
    }
    
    a <- factorial(r)*fr / sum(x[x >= r]^r)
    b <- ((n - r)*fr / ((n - r)*fr + (r + 1)*fs))^r
    rC <- 1 - a*b
    
    return(rC)
}


# Calculate first order coverage
#
# @param    x  a numeric vector of species abundance as counts
#
# @returns  Coverage estimate.
chao1Coverage <- function(x) {
    x <- x[x >= 1]
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    
    if (f2 > 0) {
        rC1 <- 1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
    } else {
        rC1 <- 1 - (f1 / n) * (((n - 1) * (f1 - 1)) / ((n - 1) * (f1 - 1) + 2))
    }
    
    return(rC1)
}


# Adjustement to observed relative abundances
#
# @param    x  vector of observed abundance counts
#
# @return   An adjusted observed species relative abundance distribution.
adjustObservedAbundance <- function(x) {
    x <- x[x >= 1]
    n <- sum(x)
    
    # Coverage
    rC1 <- calcCoverage(x, r=1)
    
    # Calculate tuning parameter
    lambda <- (1 - rC1) / sum(x/n * exp(-x))
    
    # Define adjusted relative abundance
    p <- x/n * (1 -  lambda * exp(-x))
    
    return(p)
}


# Define undetected species relative abundances
#
# @param    x  vector of detected species abundance counts.
# 
# @return   An adjusted detected species relative abundance distribution.
inferUnseenAbundance <- function(x) {
    x <- x[x >= 1]
    
    # Coverage
    rC1 <- calcCoverage(x, r=1)
    
    # Unseen count
    f0 <- inferUnseenCount(x)
    
    # Assign unseen relative abundance
    p <- rep((1 - rC1) / f0, f0)
    
    return(p)
}



#' Calculate the diversity index
#' 
#' \code{calcDiversity} calculates the clonal diversity index for a vector of diversity 
#' orders.
#'
#' @param    p  numeric vector of clone (species) counts or proportions.
#' @param    q  numeric vector of diversity orders.
#' 
#' @return   A vector of diversity scores \eqn{D} for each \eqn{q}.
#' 
#' @details
#' This method, proposed by Hill (Hill, 1973), quantifies diversity as a smooth function 
#' (\eqn{D}) of a single parameter \eqn{q}. Special cases of the generalized diversity 
#' index correspond to the most popular diversity measures in ecology: species richness 
#' (\eqn{q = 0}), the exponential of the Shannon-Weiner index (\eqn{q} approaches \eqn{1}), the 
#' inverse of the Simpson index (\eqn{q = 2}), and the reciprocal abundance of the largest 
#' clone (\eqn{q} approaches \eqn{+\infty}). At \eqn{q = 0} different clones weight equally, 
#' regardless of their size. As the parameter \eqn{q} increase from \eqn{0} to \eqn{+\infty} 
#' the diversity index (\eqn{D}) depends less on rare clones and more on common (abundant) 
#' ones, thus encompassing a range of definitions that can be visualized as a single curve. 
#' 
#' Values of \eqn{q < 0} are valid, but are generally not meaningful. The value of \eqn{D} 
#' at \eqn{q=1} is estimated by \eqn{D} at \eqn{q=0.9999}. 
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#' }
#' 
#' @seealso  Used by \code{\link{resampleDiversity}} and \code{\link{testDiversity}}.
#' 
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
#' \code{resampleDiversity} divides a set of clones by a group annotation,
#' uniformly resamples the sequences from each group, and calculates diversity
#' scores (\eqn{D}) over an interval of diversity orders (\eqn{q}).
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    min_q     minimum value of \eqn{q}.
#' @param    max_q     maximum value of \eqn{q}.
#' @param    step_q    value by which to increment \eqn{q}.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} the maximum
#'                     if automatically determined from the size of the largest group.
#' @param    replace   if \code{TRUE} resample with replacement; if \code{FALSE} resample
#'                     without replacement.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#' 
#' @return   A \code{\link{DiversityCurve}} object summarizing the diversity scores.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index proposed by 
#' Hill (Hill, 1973). See \code{\link{calcDiversity}} for further details.
#' 
#' To generate a smooth curve, \eqn{D} is calculated for each value of \eqn{q} from
#' \code{min_q} to \code{max_q} incremented by \code{step_q}; \eqn{D} at \eqn{q=1} is 
#' estimated by \eqn{D} at \eqn{q=0.9999}. Variability in total sequence counts across 
#' unique values in the \code{group} column is corrected using repeated sampling 
#' \code{nboot} times. Each resampling realization is fixed to a uniform count for each 
#' group, determined by either the value of \code{max_n} or the minimum number of 
#' sequences among all groups when \code{max_n=NULL}. 
#' 
#' The diversity index (\eqn{D}) for each group is the median value of over all resampling 
#' realizations, with the confidence interval estimate derived via the percentile method.
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#'   \item  Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
#'            peripheral blood IgE repertoires in patients with allergic rhinitis. 
#'            J Allergy Clin Immunol. 2014 134(3):604-12.
#' }
#'  
#' @seealso  See \code{\link{calcDiversity}} for the basic calculation and 
#'           \code{\link{DiversityCurve}} for the return object. 
#'           See \code{\link{testDiversity}} for significance testing.
#'           See \code{\link{plotDiversityCurve}} for plotting the return object.
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # All groups pass default minimum sampling threshold of 10 sequences
#' # With replacement
#' resampleDiversity(df, "SAMPLE", step_q=1, max_q=10, nboot=100, replace=TRUE)
#' # Without replacement
#' resampleDiversity(df, "SAMPLE", step_q=1, max_q=10, nboot=100, replace=FALSE)
#' 
#' # Increasing threshold results in exclusion of small groups and a warning message
#' resampleDiversity(df, "ISOTYPE", min_n=40, step_q=1, max_q=10, nboot=100)
#'
#' @export
resampleDiversity <- function(data, group, clone="CLONE", min_q=0, max_q=32, step_q=0.05, 
                               min_n=10, max_n=NULL, replace=FALSE, ci=0.95, nboot=2000) {
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
    pb <- txtProgressBar(min=0, max=length(group_keep), initial=0, width=40, style=3)
    div_list <- list()
    i <- 0
    for (g in group_keep) {
        i <- i + 1
        r <- which(data[[group]] == g)
        sample_mat <- replicate(nboot, data[[clone]][sample(r, n, replace=replace)])
        boot_mat <- apply(sample_mat, 2, function(x) calcDiversity(table(x), q))
        boot_ci <- t(apply(boot_mat, 1, quantile, probs=ci_probs))
        div_list[[g]] <- matrix(c(q, boot_ci[, c(2, 1, 3)]),
                                nrow=length(q), ncol=4, 
                                dimnames=list(NULL, c("q", "D", "lower", "upper")))
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    
    # Generate return object
    div <- new("DiversityCurve", 
               data=ldply(div_list, .id="group"), 
               groups=group_keep, 
               n=n, 
               nboot=nboot, 
               ci=ci)

    return(div)
}


#' Generate a clonal diversity index curve
#'
#' \code{resampleDiversity} divides a set of clones by a group annotation,
#' uniformly resamples the sequences from each group, and calculates diversity
#' scores (\eqn{D}) over an interval of diversity orders (\eqn{q}).
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    method    string defining the rarefaction method. One of \code{"depth"}, which
#'                     rarefies to the same number of sequences for each group, or 
#'                     \code{"coverage"}, which rarefies to the same degree of sample 
#'                     completeness. 
#' @param    min_q     minimum value of \eqn{q}.
#' @param    max_q     maximum value of \eqn{q}.
#' @param    step_q    value by which to increment \eqn{q}.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} the maximum
#'                     if automatically determined from the size of the largest group.
#' @param    min_c     minimum coverage required to retain a group. A group with lower coverage 
#'                     \code{min_c} prior to rarefaction will be excluded.
#' @param    replace   if \code{TRUE} resample with replacement; if \code{FALSE} resample
#'                     without replacement.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#' 
#' @return   A \code{\link{DiversityCurve}} object summarizing the diversity scores.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index proposed by 
#' Hill (Hill, 1973). See \code{\link{calcDiversity}} for further details.
#' 
#' To generate a smooth curve, \eqn{D} is calculated for each value of \eqn{q} from
#' \code{min_q} to \code{max_q} incremented by \code{step_q}.  Variability in total 
#' sequence counts across unique values in the \code{group} column is corrected by
#' rarefying to either a common number of sequences (\code{method="depth"}) or a common
#' level of sample coverage (\code{method="coverage"}) as implemented in the 
#' \code{\link[iNEXT]{iNEXT}} package and described by Chao et al, 2014. 
#' 
#' Confidence intervals are calculated using the \code{\link[iNEXT]{iNEXT}} 
#' implementation.
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45–67.
#' }
#'  
#' @seealso  See \code{\link{calcDiversity}} for the basic calculation and 
#'           \code{\link{DiversityCurve}} for the return object. 
#'           See \code{\link{testDiversity}} for significance testing.
#'           See \code{\link{plotDiversityCurve}} for plotting the return object.
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # All groups do not pass default minimum thresholds using depth rarefaction
#' div <- rarefyDiversity(df, "SAMPLE", method="depth", step_q=1, max_q=10, nboot=100)
#' plotDiversityCurve(div, legend_title="Sample")
#'                    
#' # Difference groups fail the minimum count threshold when allowing for very low coverage samples
#' div <- rarefyDiversity(df, "SAMPLE", method="coverage", step_q=1, max_q=10, min_c=0.1, nboot=100)
#' plotDiversityCurve(div, legend_title="Sample")
#'                    
#' # Grouping by isotype rather than sample identifier
#' div <- rarefyDiversity(df, "ISOTYPE", method="depth", min_n=40, step_q=1, max_q=10, nboot=100)
#' plotDiversityCurve(div, legend_title="Isotype")
#'
#' @export
rarefyDiversity <- function(data, group, clone="CLONE", method=c("depth", "coverage"), 
                            min_q=0, max_q=32, step_q=0.05, min_n=10, max_n=NULL, 
                            min_c=0.2, ci=0.95, nboot=2000) {
    #group="SAMPLE"; clone="CLONE"; method="depth"; min_q=0; max_q=32; step_q=0.05; min_n=10; max_n=NULL; min_c=0.2; ci=0.95; nboot=200
    
    # Check arguments
    method <- match.arg(method)
    
    # Verify data
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    
    if (!(group %in% names(data))) {
        stop(paste("The column", group, "does not exist in the input data.frame"))
    } else {
        data[[group]] <- as.character(data[[group]])
    }

    if (!(clone %in% names(data))) {
        stop(paste("The column", clone, "does not exist in the input data.frame"))
    } else {
        data[[clone]] <- as.character(data[[clone]])
    }

    # Calculate clonal abundance
    # TODO:  Use repertoire::countClones and implement copy number approach
    #clone_tab <- ddply(data, c(group, clone), here(summarize), COUNT=length(eval(parse(text=clone))))
    
    # Tabulate clonal abundance
    clone_tab <- data %>% 
                 group_by_(.dots=c(group, clone)) %>%
                 dplyr::summarize(COUNT=n())
    
    # Count observations per group and set sampling criteria
    #cover_tab <- ddply(clone_tab, c(group), summarize, 
    #                   COVERAGE=iNEXT:::Chat.Ind(COUNT, sum(COUNT, na.rm=TRUE)))
    #group_tab <- ddply(clone_tab, c(group), summarize, 
    #                   COUNT=sum(COUNT, na.rm=TRUE))
    #group_tab <- join(group_tab, cover_tab, by=group)

    group_tab <- clone_tab %>%
                 group_by_(.dots=c(group)) %>%
                 dplyr::summarize(SEQUENCES=sum(COUNT, na.rm=TRUE), 
                                  COVERAGE=iNEXT:::Chat.Ind(COUNT, SEQUENCES))
    group_all <- as.character(group_tab[[group]])
    group_tab <- filter(group_tab, SEQUENCES >= min_n, COVERAGE >= min_c)
    group_keep <- as.character(group_tab[[group]])
    
    if (method == "depth") {
        nsam <- min(group_tab$SEQUENCES[group_tab$SEQUENCES >= min_n], max_n)
        nsam <- setNames(rep(nsam, length(group_keep)), group_keep)
    } else if (method == "coverage") {
        # Get lowest coverage
        csam <- max(min(group_tab$COVERAGE), min_c)

        # Infer depth from coverage curves
        .fC <- function(m, x) { abs(csam - iNEXT:::Chat.Ind(x, m)) }
        nsam <- setNames(numeric(length(group_keep)), group_keep)
        for (g in group_keep) {
            x <- clone_tab$COUNT[clone_tab[[group]] == g]
            m <- group_tab$SEQUENCES[group_tab[[group]] == g]
            #y <- iNEXT:::Chat.Ind(x, 1:m)
            #nsam[g] <- which.min(abs(csam - y))
            nsam[g] <- floor(optimize(.fC, interval=c(1, m), x=x)$minimum)
        }
        # Filter to samples with enough sequences after coverage rarefaction
        nsam <- nsam[nsam >= min_n]
        group_keep <- names(nsam)
    }
    
    # Set diversity orders and confidence interval
    q <- seq(min_q, max_q, step_q)
    ci_z <- ci + (1 - ci) / 2
    
    # Warn if groups removed
    if (length(group_keep) < length(group_all)) {
        warning("Not all groups passed thresholds min_n=", min_n, " and min_c=", min_c, 
                ". Excluded: ", paste(setdiff(group_all, group_keep), collapse=", "))
    }
    
    # Generate diversity index and confidence intervals via resampling
    cat("-> CALCULATING DIVERSITY\n")
    pb <- txtProgressBar(min=0, max=length(group_keep), initial=0, width=40, style=3)
    div_list <- list()
    coverage <- setNames(numeric(length(group_keep)), group_keep)
    i <- 0
    for (g in group_keep) {
        i <- i + 1
        n <- nsam[g]
        
        # Calculate observed diversity
        abund_obs <- clone_tab$COUNT[clone_tab[[group]] == g]
        div_obs <- calcRarefiedDiversity(abund_obs, q, n) 
        #div_obs <- sapply(q, function(x) iNEXT:::Dqhat.Ind(abund_obs, q=x, m=n))
        coverage[g] <- iNEXT:::Chat.Ind(abund_obs, n)
        
        # Bootstrap diversity
        abund_inf <- iNEXT:::EstiBootComm.Ind(abund_obs)
        sample_mat <- rmultinom(nboot, n, abund_inf)
        #div_boot <- apply(sample_mat, 2, function(y) sapply(q, function(x) iNEXT:::Dqhat.Ind(y, q=x, m=n)))
        div_boot <- apply(sample_mat, 2, calcRarefiedDiversity, q=q, n=n)
        
        # Assign confidence intervals based on variance of bootstrap realizations
        sd_boot <- apply(div_boot, 1, sd)
        err_boot <- qnorm(ci_z) * sd_boot
        div_lower <- pmax(div_obs - err_boot, 0)
        div_upper <- div_obs + err_boot
        
        # Build result matrix object
        div_list[[g]] <- as.data.frame(matrix(c(q, div_obs, div_lower, div_upper),
                                       nrow=length(q), ncol=4, 
                                       dimnames=list(NULL, c("q", "D", "lower", "upper"))))
        
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    
    # TODO: Allow dplyr::tbl class for data slot of DiversityCurve
    # Generate return object
    div <- new("DiversityCurve", 
               data=as.data.frame(bind_rows(div_list, .id="group")), 
               groups=group_keep, 
               n=nsam,
               coverage=coverage,
               nboot=nboot, 
               ci=ci)
    
    return(div)
}


#' Pairwise test of the diversity index
#' 
#' \code{testDiversity} performs pairwise significance tests of the diversity index 
#' (\eqn{D}) at a given diversity order (\eqn{q}) for a set of annotation groups via
#' bootstrapping.
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
#' 
#' @return   A \code{\link{DiversityTest}} object containing p-values and summary statistics.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index proposed by 
#' Hill (Hill, 1973). See \code{\link{calcDiversity}} for further details.
#' 
#' Variability in total sequence counts across unique values in the \code{group} column is 
#' corrected using repeated subsampling with replacement (bootstrapping) \code{nboot} times. 
#' Each bootstrap realization is fixed to a uniform count for each group, determined by 
#' either the value of \code{max_n} or the minimum number of sequences among all groups when 
#' \code{max_n=NULL}. The diversity index estimate (\eqn{D}) for each group is the median value 
#' of over all bootstrap realizations. 
#' 
#' Significance of the difference in diversity index (\eqn{D}) between groups is tested by 
#' constructing a bootstrap delta distribution for each pair of unique values in the 
#' \code{group} column. The bootstrap delta distribution is built by subtracting the diversity 
#' index \eqn{Da} in \eqn{group-a} from the corresponding value \eqn{Db} in \eqn{group-b}, 
#' for all bootstrap realizations, yeilding a distribution of \code{nboot} total deltas; where 
#' \eqn{group-a} is the group with the greater median \eqn{D}. The p-value for hypothesis 
#' \eqn{Da  !=  Db} is the value of \eqn{P(0)} from the empirical cumulative distribution function of the 
#' bootstrap delta distribution, multipled by 2 for the two-tailed correction.
#' 
#' @note
#' This method may inflate statistical significance when clone sizes are uniformly small,
#' such as when most clones sizes are 1, sample size is small, and \code{max_n} is near
#' the total count of the smallest data group. Use caution when interpreting the results 
#' in such cases. We are currently investigating this potential problem.
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#'   \item  Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
#'            peripheral blood IgE repertoires in patients with allergic rhinitis. 
#'            J Allergy Clin Immunol. 2014 134(3):604-12.
#' }
#' 
#' @seealso  See \code{\link{calcDiversity}} for the basic calculation and 
#'           \code{\link{DiversityTest}} for the return object. 
#'           See \code{\link{resampleDiversity}} for curve generation.
#'           See \code{\link{ecdf}} for computation of the empirical cumulative 
#'           distribution function.
#' 
#' @examples  
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Groups under the size threshold are excluded and a warning message is issued.
#' testDiversity(df, "SAMPLE", q=0, min_n=30, nboot=1000)
#' 
#' @export
testDiversity <- function(data, q, group, clone="CLONE", min_n=10, max_n=NULL, nboot=2000) {
    
    # TODO:  write plotDiversityTest function
    
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
    pb <- txtProgressBar(min=0, max=ngroup, initial=0, width=40, style=3)
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
    group_pairs <- combn(group_keep, 2, simplify=F)
    npairs <- length(group_pairs)
    delta_mat <- matrix(NA, nboot, npairs)
    pvalue_mat <- matrix(NA, npairs, 5, dimnames=list(NULL, c("pvalue", 
                                                              "delta_median", 
                                                              "delta_mad", 
                                                              "delta_mean", 
                                                              "delta_sd")))
    test_names <- sapply(group_pairs, paste, collapse=" != ")
    for (i in 1:npairs) {
        g1 <- group_pairs[[i]][1]
        g2 <- group_pairs[[i]][2]
        # TODO:  verify this is correct. Is g1 - g2 different from g2 - g1?
        if (median(div_mat[, g1]) >= median(div_mat[, g2])) {
            g_delta <- div_mat[, g1] - div_mat[, g2]
        } else {
            g_delta <- div_mat[, g2] - div_mat[, g1]
        }  
        
        # Determine p-value
        g_cdf <- ecdf(g_delta)
        p <- g_cdf(0)
        p <- ifelse(p <= 0.5, p * 2, (1 - p) * 2)
        pvalue_mat[i, ] <- c(p, 
                             median(g_delta), 
                             mad(g_delta), 
                             mean(g_delta), 
                             sd(g_delta))
    }
    
    tests_df <- cbind(data.frame(test=test_names), as.data.frame(pvalue_mat))
    summary_df <- data.frame(group=group_keep, 
                             median=apply(div_mat, 2, median),
                             mad=apply(div_mat, 2, mad),
                             mean=apply(div_mat, 2, mean),
                             sd=apply(div_mat, 2, sd))
    
    # Generate return object
    div <- new("DiversityTest", 
               tests=tests_df, 
               summary=summary_df,
               groups=group_keep,
               q=q,
               n=n, 
               nboot=nboot)
    
    return(div)
}


#### Plotting functions ####

#' Plot the results of resampleDiversity
#' 
#' \code{plotDiversityCurve} plots a \code{DiversityCurve} object.
#'
#' @param    data            \code{\link{DiversityCurve}} object returned by 
#'                           \code{\link{rarefyDiversity}}.
#' @param    colors          named character vector whose names are values in the 
#'                           \code{group} column of the \code{data} slot of \code{data},
#'                           and whose values are colors to assign to those group names.
#' @param    main_title      string specifying the plot title.
#' @param    legend_title    string specifying the legend title.
#' @param    log_q           if \code{TRUE} then plot \eqn{q} on a log scale;
#'                           if \code{FALSE} plot on a linear scale.
#' @param    log_d           if \code{TRUE} then plot the diversity scores \eqn{D} 
#'                           on a log scale; if \code{FALSE} plot on a linear scale.
#' @param    xlim            numeric vector of two values specifying the 
#'                           \code{c(lower, upper)} x-axis limits.
#' @param    ylim            numeric vector of two values specifying the 
#'                           \code{c(lower, upper)} y-axis limits.
#' @param    annotate        string defining whether to added values to the group labels 
#'                           of the legend showing sequence counts (\code{"depth"}),
#'                           coverage (\code{"coverage"}), both counts and coverage 
#'                           (\code{"both"}), or neither (\code{"none"}).
#' @param    silent          if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                           object; if \code{FALSE} draw the plot.
#' @param    ...             additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  See \code{\link{rarefyDiversity}} for generating \code{\link{DiversityCurve}}
#'           objects for input. Plotting is performed with \code{\link{ggplot}}.
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # All groups pass default minimum sampling threshold of 10 sequences
#' div <- rarefyDiversity(df, "SAMPLE", min_c=0.1, step_q=0.1, max_q=10, nboot=100)
#' plotDiversityCurve(div, legend_title="Sample")
#' 
#' @export
plotDiversityCurve <- function(data, colors=NULL, main_title="Diversity", 
                               legend_title=NULL, log_q=TRUE, log_d=TRUE,
                               xlim=NULL, ylim=NULL, 
                               annotate=c("both", "depth", "coverage", "none"),
                               silent=FALSE, ...) {
    # Check arguments
    annotate <- match.arg(annotate)
    
    # Define group label annotations
    if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        group_labels <- setNames(paste0(data@groups, 
                                        " (N=", data@n, ")"), 
                                 data@groups)
    } else if (annotate == "coverage") {
        group_labels <- setNames(paste0(data@groups, 
                                        " (C=", round(data@coverage, 3), ")"), 
                                 data@groups)
    } else if (annotate == "both") {
        group_labels <- setNames(paste0(data@groups, 
                                        " (N=", data@n, 
                                        ", C=", round(data@coverage, 3), ")"), 
                                 data@groups)
    }
    
    # Define base plot elements
    p1 <- ggplot(data@data, aes(x=q, y=D, group=group)) + 
        ggtitle(main_title) + 
        getBaseTheme() + 
        xlab('q') +
        ylab(expression(''^q * D)) +
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=0.4) +
        geom_line(aes(color=group))
    
    # Set colors and legend
    if (!is.null(colors)) {
        p1 <- p1 + scale_color_manual(name=legend_title, labels=group_labels, values=colors) +
            scale_fill_manual(name=legend_title, labels=group_labels, values=colors)
    } else {
        p1 <- p1 + scale_color_discrete(name=legend_title, labels=group_labels) +
            scale_fill_discrete(name=legend_title, labels=group_labels)
    }
    
    # Set x-axis style
    if (log_q) {
        p1 <- p1 + scale_x_continuous(trans=log2_trans(), limits=xlim,
                                      breaks=trans_breaks('log2', function(x) 2^x),
                                      labels=trans_format('log2', math_format(2^.x)))
    } else {
        p1 <- p1 + scale_x_continuous(limits=xlim)
    }
    
    # Set y-axis style
    if (log_d) {
        p1 <- p1 + scale_y_continuous(trans=log2_trans(), limits=ylim,
                                      breaks=trans_breaks('log2', function(x) 2^x),
                                      labels=trans_format('log2', math_format(2^.x)))
    } else {
        p1 <- p1 + scale_y_continuous(limits=ylim)
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))

    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}