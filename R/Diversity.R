# Clonal diversity analysis

#' @include Classes.R
NULL


#### Coverage functions ####

#' Calculate sample coverage
#' 
#' \code{calcCoverage} calculates the sample coverage estimate, a measure of sample 
#' completeness, for varying orders using the method of Chao et al, 2015, falling back 
#' to the Chao1 method in the first order case.
#'
#' @param    x  numeric vector of abundance counts.
#' @param    r  coverage order to calculate.
#' 
#' @return   The sample coverage of the given order \code{r}.
#' 
#' @references
#' \enumerate{
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#' 
#' @seealso  
#' Used by \link{rarefyDiversity}.
#'           
#' @examples
#' # Calculate clone sizes
#' clones <- countClones(ExampleDb, groups="SAMPLE")
#' 
#' # Calculate 1first order coverage for a single sample
#' calcCoverage(clones$SEQ_COUNT[clones$SAMPLE == "+7d"])
#'
#' @export
calcCoverage <- function(x, r=1) {
    # Use traditional calculation for 1st order coverage
    if (r == 1) { return(calcChao1Coverage(x)) }
    
    # Use general form for 2nd order and higher coverage
    x <- x[x >= 1]
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
calcChao1Coverage <- function(x) {
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


# Calculates diversity under rarefaction
# 
# Calculates Hill numbers under rarefaction
#
# @param    x  vector of observed abundance counts.
# @param    m  the sequence count to rarefy to.
#
# @return   The first order coverage estimate
inferRarefiedCoverage <- function(x, m) {
    x <- x[x >= 1]
    n <- sum(x)
    if (m > n) {
        stop("m must be <= the total count of observed sequences.")
    }
    
    # Unrarefied case
    if (m == n) {
        return(calcCoverage(x, r=1))
    }
    
    # Calculate rarefied coverage
    # TODO: Read up on this and fix
    #rC1 <- iNEXT:::Chat.Ind(x, m)
    y <- x[(n - x) >= m]
    rC1 <- 1 - sum(y/n * exp(lgamma(n - y + 1) - lgamma(n - m - y + 1) - lgamma(n) + lgamma(n - m)))
    
    return(rC1)
}


#### Abundance functions ####

# Calculate undetected species
# 
# Calculates the lower bound of undetected species counts using the Chao1 estimator.
#
# @param    x  vector of observed abundance counts.
# 
# @return   The count of undetected species.
inferUnseenCount <- function(x) {
    x <- x[x >= 1]
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


#' Tabulates clones sizes
#' 
#' \code{countClones} determines the number of sequences and total copy number of 
#' clonal groups.
#'
#' @param    data    data.frame with Change-O style columns containing clonal assignments.
#' @param    groups  character vector defining \code{data} columns containing grouping 
#'                   variables. If \code{group=NULL}, then do not group data.
#' @param    copy    name of the \code{data} column containing copy numbers for each 
#'                   sequence. If this value is specified, then total copy abundance
#'                   is determined by the sum of copy numbers within each clonal group.
#' @param    clone   name of the \code{data} column containing clone identifiers.
#' 
#' @return   A data.frame summarizing clone counts and frequencies with columns:
#'           \itemize{
#'             \item \code{CLONE}:       clone identifier.
#'             \item \code{SEQ_COUNT}:   total number of sequences for the clone.
#'             \item \code{SEQ_FREQ}:    frequency of the clone as a fraction of the total
#'                                       number of sequences within each group.
#'             \item \code{COPY_COUNT}:  sum of the copy counts in the \code{copy} column.
#'                                       Only present if the \code{copy} argument is 
#'                                       specified.
#'             \item \code{COPY_FREQ}:   frequency of the clone as a fraction of the total
#'                                       copy number within each group. Only present if 
#'                                       the \code{copy} argument is specified.
#'           }
#'           Also includes additional columns specified in the \code{groups} argument.
#' 
#' @examples
#' # Without copy numbers
#' clones <- countClones(ExampleDb, groups="SAMPLE")
#'
#' # With copy numbers and multiple groups
#' clones <- countClones(ExampleDb, groups=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")
#' 
#' @export
countClones <- function(data, groups=NULL, copy=NULL, clone="CLONE") {
    # Check input
    check <- checkColumns(data, c(clone, copy, groups))
    if (check != TRUE) { stop(check) }
    
    # Tabulate clonal abundance
    if (is.null(copy)) {
        clone_tab <- data %>% 
            group_by_(.dots=c(groups, clone)) %>%
            dplyr::summarize(SEQ_COUNT=n()) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(SEQ_COUNT)") %>%
            dplyr::rename_(.dots=c("CLONE"=clone))
    } else {
        clone_tab <- data %>% 
            group_by_(.dots=c(groups, clone)) %>%
            dplyr::summarize_(SEQ_COUNT=interp(~length(x), x=as.name(clone)),
                              COPY_COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy))) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT")),
                           COPY_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("COPY_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(COPY_COUNT)") %>%
            dplyr::rename_(.dots=c("CLONE"=clone))
    }
    
    return(clone_tab)
}

# Perform boostrap abundance calculation
# 
# @param    x      named vector of observed abundance values.
# @param    n      number of samples to draw from the estimate complete abundance distribution.
# @param    z      Confidence interval tail probability. Meaning, \code{ci + (1 - ci) / 2)}.
# @param    nboot  number of bootstrap realizations.
# 
# @return   A data.frame summarizing abundance with columns:
#           \itemize{
#             \item \code{CLONE}:  clone identifier.
#             \item \code{P}:      abundance (clonal frequency).
#             \item \code{LOWER}:  lower confidence interval bound.
#             \item \code{UPPER}:  upper confidence interval bound.
#             \item \code{RANK}:   adundance rank. Lower ranks correspond to higher abundance values.
#                                  Ties of P get increasing rank values.
#                                  
bootstrapAbundance <- function(x, n, z=0.975, nboot=2000) {
    # Calculate estimated complete abundance distribution
    p1 <- adjustObservedAbundance(x)
    p2 <- inferUnseenAbundance(x)
    p <- c(p1, p2)
    p2_names <- if (length(p2) > 0) { paste0("U", 1:length(p2)) } else { NULL }
    names(p) <- c(names(x), p2_names)
    
    # Bootstrap abundance
    boot_mat <- rmultinom(nboot, n, p) / n
    
    # Assign confidence intervals based on variance of bootstrap realizations
    boot_sd <- apply(boot_mat, 1, sd)
    boot_err <- qnorm(z) * boot_sd
    p_lower <- pmax(p - boot_err, 0)
    p_upper <- p + boot_err
    
    # Assemble and sort abundance data.frame
    abund_df <- tibble(CLONE=names(p), P=p, LOWER=p_lower, UPPER=p_upper)
    abund_df <- dplyr::arrange_(abund_df, .dots="desc(P)")
    abund_df$RANK <- 1:nrow(abund_df)
    
    return(abund_df)
}


#' Estimates the complete clonal relative abundance distribution
#' 
#' \code{estimateAbundance} estimates the complete clonal relative abundance distribution 
#' and confidence intervals on clone sizes using bootstrapping.
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    group     name of the \code{data} column containing group identifiers. 
#'                     If \code{NULL} then no grouping is performed and the \code{GROUP} 
#'                     column of the output will contain the value \code{NA} for each row.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each 
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of 
#'                     copy numbers within each clonal group.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#' @param    progress  if \code{TRUE} show a progress bar. 
#' 
#' @return   A data.frame with relative clonal abundance data and confidence intervals,
#'           containing the following columns:
#'           \itemize{
#'             \item  \code{GROUP}:  group identifier. Will be code{NA} if \code{group=NULL}.
#'             \item  \code{CLONE}:  clone identifier.
#'             \item  \code{P}:      relative abundance of the clone.
#'             \item  \code{LOWER}:  lower confidence inverval bound.
#'             \item  \code{UPPER}:  upper confidence interval bound.
#'             \item  \code{RANK}:   the rank of the clone abundance.
#'           }
#'           
#' @details
#' The complete clonal abundance distribution determined inferred by using the Chao1 
#' estimator to estimate the number of seen clones, and then applying the relative abundance 
#' correction and unseen clone frequencies described in Chao et al, 2015.
#'
#' Confidence intervals are derived using the standard deviation of the resampling 
#' realizations, as described in Chao et al, 2015.
#' 
#' @references
#' \enumerate{
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#' 
#' @seealso  
#' See \link{plotAbundanceCurve} for plotting of the abundance distribution.
#' See \link{rarefyDiversity} for a similar application to clonal diversity.
#'           
#' @examples
#' abund <- estimateAbundance(ExampleDb, "SAMPLE", nboot=100)
#'
#' @export
estimateAbundance <- function(data, group=NULL, clone="CLONE", copy=NULL, ci=0.95, 
                              nboot=2000, progress=FALSE) {
    ## DEBUG
    # group="SAMPLE"; clone="CLONE"; copy="UID_CLUSTCOUNT"; ci=0.95; nboot=200
    # data=clones; group="SUBJECT"; clone="CLONE"; copy=NULL; ci=0.95; nboot=200; progress=FALSE
    
    # Check input
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    check <- checkColumns(data, c(group, clone, copy))
    if (check != TRUE) { stop(check) }
    
    # Tabulate clonal abundance
    if (is.null(copy)) {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize(COUNT=n())
    } else {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize_(COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy)))
    }
    
    # Set confidence interval
    ci_z <- ci + (1 - ci) / 2
    
    if (!is.null(group)) {
        # Summarize groups
        group_tab <- clone_tab %>%
            group_by_(.dots=c(group)) %>%
            dplyr::summarize_(SEQUENCES=interp(~sum(x, na.rm=TRUE), x=as.name("COUNT")))
        group_set <- as.character(group_tab[[group]])
        nsam <- setNames(group_tab$SEQUENCES, group_set)

        # Generate diversity index and confidence intervals via resampling
        if (progress) { 
            pb <- progressBar(length(group_set))
        }
        abund_list <- list()
        for (g in group_set) {
            n <- nsam[g]
            
            # Extract abundance vector
            abund_obs <- clone_tab$COUNT[clone_tab[[group]] == g]
            names(abund_obs) <- clone_tab[[clone]][clone_tab[[group]] == g]
     
            # Infer complete abundance distribution
            abund_list[[g]] <- bootstrapAbundance(abund_obs, n, z=ci_z, nboot=nboot)
            
            if (progress) { pb$tick() }
        }
        curve_df <- as.data.frame(bind_rows(abund_list, .id="GROUP"))
    } else {
        # Summarize counts
        group_set <- as.character(NA)
        nsam <- sum(clone_tab$COUNT, na.rm=TRUE)
        
        # Extract abundance vector
        abund_obs <- clone_tab$COUNT
        names(abund_obs) <- clone_tab[[clone]]
        
        # Infer complete abundance distribution
        curve_df <- bootstrapAbundance(abund_obs, nsam, z=ci_z, nboot=nboot)
        curve_df$GROUP <- NA
    }
    # Generate return object
    curve <- new("AbundanceCurve", 
                 data=curve_df, 
                 groups=group_set, 
                 n=nsam,
                 nboot=nboot, 
                 ci=ci)
    
    return(curve)
}


#### Diversity functions ####

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
#' @seealso  Used by \link{rarefyDiversity} and \link{testDiversity}.
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
    # Add jitter to q=1
    q[q == 1] <- 0.9999
    # Remove zeros
    p <- p[p > 0]
    # Convert p to proportional abundance
    p <- p / sum(p)
    # Calculate D for each q
    D <- sapply(q, function(x) sum(p^x)^(1 / (1 - x)))
    
    return(D)
}


# Calculates diversity under rarefaction
# 
# Calculates Hill numbers under rarefaction
#
# @param    x  vector of observed abundance counts.
# @param    q  numeric vector of diversity orders.
# @param    m  the sequence count to rarefy to.
#
# @return   A vector of diversity scores \eqn{D} for each \eqn{q}.
inferRarefiedDiversity <- function(x, q, m) {
    x <- x[x >= 1]
    n <- sum(x)
    if (m > n) {
        stop("m must be <= the total count of observed sequences.")
    }
    q[q == 1] <- 0.9999
    
    # Tabulate frequency counts from 1:n
    fk_n <- tabulate(x, nbins=n)
    
    # Calculate estimated fk(m)
    fk_m <- sapply(1:m, function(k) sum(exp(lchoose(k:m, k) + 
                                            lchoose(n - k:m, m - k) - 
                                            lchoose(n, m))*fk_n[k:m]))
    
    # Calculate diversity
    D <- sapply(q, function(r) sum((1:m / m)^r * fk_m)^(1 / (1 - r)))
    
    return(D)
}


#' Generate a clonal diversity index curve
#'
#' \code{rarefyDiversity} divides a set of clones by a group annotation,
#' resamples the sequences from each group, and calculates diversity
#' scores (\eqn{D}) over an interval of diversity orders (\eqn{q}).
#' 
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each 
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of 
#'                     copy numbers within each clonal group.
#' @param    min_q     minimum value of \eqn{q}.
#' @param    max_q     maximum value of \eqn{q}.
#' @param    step_q    value by which to increment \eqn{q}.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} then no 
#'                     maximum is set.
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#' @param    uniform   if \code{TRUE} then uniformly resample each group to the same 
#'                     number of observations. If \code{FALSE} then allow each group to
#'                     be resampled to its original size or, if specified, \code{max_size}.
#' @param    progress  if \code{TRUE} show a progress bar.
#' 
#' @return   A \link{DiversityCurve} object summarizing the diversity scores.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index (Hill numbers) 
#' proposed by Hill (Hill, 1973). See \link{calcDiversity} for further details.
#'
#' Diversity is calculated on the estimated complete clonal abundance distribution.
#' This distribution is inferred by using the Chao1 estimator to estimate the number
#' of seen clones, and applying the relative abundance correction and unseen clone
#' frequency described in Chao et al, 2015.
#'
#' To generate a smooth curve, \eqn{D} is calculated for each value of \eqn{q} from
#' \code{min_q} to \code{max_q} incremented by \code{step_q}.  When \code{uniform=TRUE}
#' variability in total sequence counts across unique values in the \code{group} column 
#' is corrected by repeated resampling from the estimated complete clonal distribution to a 
#' common number of sequences.
#' 
#' The diversity index (\eqn{D}) for each group is the mean value of over all resampling 
#' realizations. Confidence intervals are derived using the standard deviation of the 
#' resampling realizations, as described in Chao et al, 2015.
#' 
#' @references
#' \enumerate{
#'   \item  Hill M. Diversity and evenness: a unifying notation and its consequences. 
#'            Ecology. 1973 54(2):427-32.
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#'  
#' @seealso  See \link{calcDiversity} for the basic calculation and 
#'           \link{DiversityCurve} for the return object. 
#'           See \link{testDiversity} for significance testing.
#'           See \link{plotDiversityCurve} for plotting the return object.
#' 
#' @examples
#' # Group by sample identifier
#' div <- rarefyDiversity(ExampleDb, "SAMPLE", step_q=1, max_q=10, nboot=100)
#' plotDiversityCurve(div, legend_title="Sample")
#'                    
#' # Grouping by isotype rather than sample identifier
#' div <- rarefyDiversity(ExampleDb, "ISOTYPE", min_n=40, step_q=1, max_q=10, 
#'                        nboot=100)
#' plotDiversityCurve(div, legend_title="Isotype")
#'
#' @export
rarefyDiversity <- function(data, group, clone="CLONE", copy=NULL, 
                            min_q=0, max_q=4, step_q=0.05, min_n=30, max_n=NULL, 
                            ci=0.95, nboot=2000, uniform=TRUE, progress=FALSE) {
    #group="SAMPLE"; clone="CLONE"; copy=NULL; min_q=0; max_q=4; step_q=1; min_n=30; max_n=NULL; ci=0.95; nboot=200

    # Check input
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    check <- checkColumns(data, c(group, clone, copy))
    if (check != TRUE) { stop(check) }
    
    # Cast grouping to columns to character
    data[[group]] <- as.character(data[[group]])
    data[[clone]] <- as.character(data[[clone]])

    # Tabulate clonal abundance
    # TODO: Can this be replaced by countClones?
    if (is.null(copy)) {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize(COUNT=n())
    } else {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize_(COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy)))
    }
    
    # Count observations per group and set sampling criteria
    group_tab <- clone_tab %>%
        group_by_(.dots=c(group)) %>%
        dplyr::summarize_(SEQUENCES=interp(~sum(x, na.rm=TRUE), x=as.name("COUNT")))
    group_all <- as.character(group_tab[[group]])
    group_tab <- group_tab[group_tab$SEQUENCES >= min_n, ]
    group_keep <- as.character(group_tab[[group]])
    
    # Set number of samples sequence
    if (uniform) {
        nsam <- min(group_tab$SEQUENCES, max_n)
        nsam <- setNames(rep(nsam, length(group_keep)), group_keep)
    } else {
        nsam <- if (is.null(max_n)) { group_tab$SEQUENCES } else { pmin(group_tab$SEQUENCES, max_n) }
        nsam <- setNames(nsam, group_keep)
    }

    # Set diversity orders and confidence interval
    q <- seq(min_q, max_q, step_q)
    ci_z <- ci + (1 - ci) / 2
    
    # Check for q=0 and set index of q=0 
    q0 <- (0 %in% q)
    if (!q0) { q <- c(0, q) }
    qi <- which(q == 0)
    
    # Warn if groups removed
    if (length(group_keep) < length(group_all)) {
        warning("Not all groups passed threshold min_n=", min_n, ".", 
                "Excluded: ", paste(setdiff(group_all, group_keep), collapse=", "))
    }
    
    # Generate diversity index and confidence intervals via resampling
    if (progress) { 
        pb <- progressBar(length(group_keep))
    }
    div_list <- list()
    #coverage <- setNames(numeric(length(group_keep)), group_keep)
    for (g in group_keep) {
        n <- nsam[g]
        
        # Get observed abundance
        abund_obs <- clone_tab$COUNT[clone_tab[[group]] == g]

        # Calculate rarefied coverage
        #coverage[g] <- iNEXT:::Chat.Ind(abund_obs, n)
        #coverage[g] <- inferRarefiedCoverage(abund_obs, n)
        
        # Estimate complete abundance distribution
        #abund_inf <- iNEXT:::EstiBootComm.Ind(abund_obs)
        p1 <- adjustObservedAbundance(abund_obs)
        p2 <- inferUnseenAbundance(abund_obs)
        abund_inf <- c(p1, p2)

        # Generate bootstrap distributions
        sample_mat <- rmultinom(nboot, n, abund_inf)
        
        # Calculate diversity and coverage
        div_boot <- apply(sample_mat, 2, calcDiversity, q=q)
        #cover_boot <- apply(sample_mat, 2, calcCoverage, r=1)
        
        # Assign observed diversity and standard deviation from bootstrap realizations
        div_obs <- apply(div_boot, 1, mean)
        div_sd <- apply(div_boot, 1, sd)
        # Diversity confidence intervals
        div_error <- qnorm(ci_z) * div_sd
        div_lower <- pmax(div_obs - div_error, 0)
        div_upper <- div_obs + div_error
        # Evenness
        even_obs <- div_obs / div_obs[qi]
        even_lower <- div_lower / div_obs[qi]
        even_upper <- div_upper / div_obs[qi]
        
        # Build result matrix        
        result_mat <- matrix(c(q, div_obs, div_sd, div_lower, div_upper,
                               even_obs, even_lower, even_upper),
                             nrow=length(q), ncol=8, 
                             dimnames=list(NULL, c("Q", "D", "D_SD", "D_LOWER", "D_UPPER",
                                                   "E", "E_LOWER", "E_UPPER")))
        
        # Remove q=0 if required
        if (!q0) { result_mat <- result_mat[-qi, ] }
        
        # Update list with results
        div_list[[g]] <- as.data.frame(result_mat)
        
        if (progress) { pb$tick() }
    }
    
    # TODO: Allow dplyr::tbl class for data slot of DiversityCurve
    # Generate return object
    div <- new("DiversityCurve", 
               data=as.data.frame(bind_rows(div_list, .id="GROUP")), 
               groups=group_keep, 
               n=nsam,
               nboot=nboot, 
               ci=ci)
    
    return(div)
}


#' Pairwise test of the diversity index
#' 
#' \code{testDiversity} performs pairwise significance tests of the diversity index 
#' (\eqn{D}) at a given diversity order (\eqn{q}) for a set of annotation groups via
#' rarefaction and bootstrapping.
#'
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    q         diversity order to test.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each 
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of 
#'                     copy numbers within each clonal group.
#' @param    min_n     minimum number of observations to sample.
#'                     A group with less observations than the minimum is excluded.
#' @param    max_n     maximum number of observations to sample. If \code{NULL} the maximum
#'                     if automatically determined from the size of the largest group.
#' @param    nboot     number of bootstrap realizations to perform.
#' @param    progress  if \code{TRUE} show a progress bar.
#' 
#' @return   A \link{DiversityTest} object containing p-values and summary statistics.
#' 
#' @details
#' Clonal diversity is calculated using the generalized diversity index proposed by 
#' Hill (Hill, 1973). See \link{calcDiversity} for further details.
#' 
#' Diversity is calculated on the estimated complete clonal abundance distribution.
#' This distribution is inferred by using the Chao1 estimator to estimate the number
#' of seen clones, and applying the relative abundance correction and unseen clone
#' frequency described in Chao et al, 2014.
#'
#' Variability in total sequence counts across unique values in the \code{group} column is 
#' corrected by repeated resampling from the estimated complete clonal distribution to 
#' a common number of sequences. The diversity index estimate (\eqn{D}) for each group is 
#' the mean value of over all bootstrap realizations. 
#' 
#' Significance of the difference in diversity index (\eqn{D}) between groups is tested by 
#' constructing a bootstrap delta distribution for each pair of unique values in the 
#' \code{group} column. The bootstrap delta distribution is built by subtracting the diversity 
#' index \eqn{Da} in \eqn{group-a} from the corresponding value \eqn{Db} in \eqn{group-b}, 
#' for all bootstrap realizations, yeilding a distribution of \code{nboot} total deltas; where 
#' \eqn{group-a} is the group with the greater mean \eqn{D}. The p-value for hypothesis 
#' \eqn{Da  !=  Db} is the value of \eqn{P(0)} from the empirical cumulative distribution 
#' function of the bootstrap delta distribution, multiplied by 2 for the two-tailed correction.
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
#'   \item  Chao A. Nonparametric Estimation of the Number of Classes in a Population. 
#'            Scand J Stat. 1984 11, 265270.
#'   \item  Wu Y-CB, et al. Influence of seasonal exposure to grass pollen on local and 
#'            peripheral blood IgE repertoires in patients with allergic rhinitis. 
#'            J Allergy Clin Immunol. 2014 134(3):604-12.
#'   \item  Chao A, et al. Rarefaction and extrapolation with Hill numbers: 
#'            A framework for sampling and estimation in species diversity studies. 
#'            Ecol Monogr. 2014 84:45-67.
#'   \item  Chao A, et al. Unveiling the species-rank abundance distribution by 
#'            generalizing the Good-Turing sample coverage theory. 
#'            Ecology. 2015 96, 11891201.
#' }
#' 
#' @seealso  See \link{calcDiversity} for the basic calculation and 
#'           \link{DiversityTest} for the return object. 
#'           See \link{rarefyDiversity} for curve generation.
#'           See \link{ecdf} for computation of the empirical cumulative 
#'           distribution function.
#' 
#' @examples  
#' # Groups under the size threshold are excluded and a warning message is issued.
#' testDiversity(ExampleDb, "SAMPLE", q=0, min_n=30, nboot=100)
#' 
#' @export
testDiversity <- function(data, q, group, clone="CLONE", copy=NULL, 
                          min_n=30, max_n=NULL, nboot=2000, progress=FALSE) {
    #group="SAMPLE"; clone="CLONE"; copy=NULL; q=1; min_n=30; max_n=NULL; nboot=200
    # TODO:  write plotDiversityTest function

    # Check input
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    check <- checkColumns(data, c(group, clone, copy))
    if (check != TRUE) { stop(check) }
    
    # Cast grouping to columns to character
    data[[group]] <- as.character(data[[group]])
    data[[clone]] <- as.character(data[[clone]])
    
    # Tabulate clonal abundance
    if (is.null(copy)) {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize(COUNT=n())
    } else {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize_(COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy)))
    }
    
    # Count observations per group and set sampling criteria
    group_tab <- clone_tab %>%
        group_by_(.dots=c(group)) %>%
        dplyr::summarize_(SEQUENCES=interp(~sum(x, na.rm=TRUE), x=as.name("COUNT")))
    group_all <- as.character(group_tab[[group]])
    group_tab <- group_tab[group_tab$SEQUENCES >= min_n, ]
    group_keep <- as.character(group_tab[[group]])
    
    # Set number of samples per group
    nsam <- min(group_tab$SEQUENCES, max_n)
    nsam <- setNames(rep(nsam, length(group_keep)), group_keep)
    
    # Warn if groups removed
    if (length(group_keep) < length(group_all)) {
        warning("Not all groups passed threshold min_n=", min_n, ". ", 
                "Excluded: ", paste(setdiff(group_all, group_keep), collapse=", "))
    }

    # Generate diversity index and confidence intervals via resampling
    if (progress) { 
        pb <- progressBar(length(group_keep))
    }
    ngroup <- length(group_keep)
    div_mat <- matrix(NA, nboot, ngroup, dimnames=list(NULL, group_keep))
    for (i in 1:ngroup) {
        g <- group_keep[i]
        m <- nsam[g]
        
        # Estimate complete abundance distribution
        #abund_inf <- iNEXT:::EstiBootComm.Ind(abund_obs)
        abund_obs <- clone_tab$COUNT[clone_tab[[group]] == g]
        p1 <- adjustObservedAbundance(abund_obs)
        p2 <- inferUnseenAbundance(abund_obs)
        abund_inf <- c(p1, p2)

        # Generate bootstrap distributions
        sample_mat <- rmultinom(nboot, m, abund_inf)
        
        # Calculate diversity
        div_mat[, i] <- apply(sample_mat, 2, calcDiversity, q=q)

        if (progress) { pb$tick() }
    }
        
    # Compute ECDF of bootstrap distribution shift from bootstrap deltas
    group_pairs <- combn(group_keep, 2, simplify=F)
    npairs <- length(group_pairs)
    pvalue_mat <- matrix(NA, npairs, 3, 
                         dimnames=list(NULL, c("DELTA_MEAN", "DELTA_SD", "PVALUE")))
    test_names <- sapply(group_pairs, paste, collapse=" != ")
    for (i in 1:npairs) {
        g1 <- group_pairs[[i]][1]
        g2 <- group_pairs[[i]][2]
        # TODO:  verify this is correct. Is g1 - g2 different from g2 - g1?
        if (mean(div_mat[, g1]) >= mean(div_mat[, g2])) {
            g_delta <- div_mat[, g1] - div_mat[, g2]
        } else {
            g_delta <- div_mat[, g2] - div_mat[, g1]
        }  
        
        # Determine p-value
        g_cdf <- ecdf(g_delta)
        p <- g_cdf(0)
        p <- ifelse(p <= 0.5, p * 2, (1 - p) * 2)
        pvalue_mat[i, ] <- c(mean(g_delta), sd(g_delta), p)
    }
    
    tests_df <- cbind(data.frame(test=test_names), as.data.frame(pvalue_mat))
    summary_df <- data.frame(GROUP=group_keep, 
                             MEAN=apply(div_mat, 2, mean),
                             SD=apply(div_mat, 2, sd))
    
    # Generate return object
    div <- new("DiversityTest", 
               tests=as.data.frame(tests_df), 
               summary=as.data.frame(summary_df),
               groups=group_keep,
               q=q,
               n=nsam, 
               nboot=nboot)
    
    return(div)
}


#### Plotting functions ####

#' Plots a clonal abundance distribution
#' 
#' \code{plotAbundanceCurve} plots the results from estimating the complete clonal 
#' relative abundance distribution. The distribution is plotted as a log rank abundance 
#' distribution.
#' 
#' @param    data          \link{AbundanceCurve} object returned by \link{estimateAbundance}.
#' @param    colors        named character vector whose names are values in the 
#'                         \code{group} column of \code{data} and whose values are 
#'                         colors to assign to those group names.
#' @param    main_title    string specifying the plot title.
#' @param    legend_title  string specifying the legend title.
#' @param    xlim          numeric vector of two values specifying the 
#'                         \code{c(lower, upper)} x-axis limits.
#' @param    ylim          numeric vector of two values specifying the 
#'                         \code{c(lower, upper)} y-axis limits.
#' @param    annotate      string defining whether to added values to the group labels 
#'                         of the legend. When \code{"none"} (default) is specified no
#'                         annotations are added. Specifying (\code{"depth"}) adds 
#'                         sequence counts to the labels.
#' @param    silent        if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                         object; if \code{FALSE} draw the plot.
#' @param    ...           additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  
#' See \link{AbundanceCurve} for the input object and \link{estimateAbundance} for 
#' generating the input abundance distribution.
#' Plotting is performed with \link{ggplot}.
#'           
#' @examples
#' # Estimate abundance by sample and plot
#' abund <- estimateAbundance(ExampleDb, "SAMPLE", nboot=100)
#' plotAbundanceCurve(abund, legend_title="Sample")
#' 
#' @export
plotAbundanceCurve <- function(data, colors=NULL, main_title="Rank Abundance", 
                               legend_title=NULL, xlim=NULL, ylim=NULL, 
                               annotate=c("none", "depth"),
                               silent=FALSE, ...) {
    # Check arguments
    annotate <- match.arg(annotate)
    
    # Define group label annotations
    if (all(is.na(data@groups))) {
        group_labels <- NA  
    } else if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        group_labels <- setNames(paste0(data@groups, " (N=", data@n, ")"), 
                                 data@groups)
    }
    
    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL
    
    if (any(!is.na(data@groups))) {
        # Define grouped plot
        p1 <- ggplot(data@data, aes_string(x="RANK", y="P", group="GROUP")) + 
            ggtitle(main_title) + 
            baseTheme() + 
            xlab('Rank') +
            ylab('Abundance') +
            scale_x_log10(limits=xlim,
                          breaks=scales::trans_breaks('log10', function(x) 10^x),
                          labels=scales::trans_format('log10', scales::math_format(10^.x))) +
            scale_y_continuous(labels=scales::percent) +
            geom_ribbon(aes_string(ymin="LOWER", ymax="UPPER", fill="GROUP"), alpha=0.4) +
            geom_line(aes_string(color="GROUP"))
        
        # Set colors and legend
        if (!is.null(colors)) {
            p1 <- p1 + scale_color_manual(name=legend_title, labels=group_labels, values=colors) +
                scale_fill_manual(name=legend_title, labels=group_labels, values=colors)
        } else {
            p1 <- p1 + scale_color_discrete(name=legend_title, labels=group_labels) +
                scale_fill_discrete(name=legend_title, labels=group_labels)
        }
    } else {
        # Set color
        if (!is.null(colors) & length(colors) == 1) {
            line_color <- colors
        } else {
            line_color <- "black"
        }
        # Define plot
        p1 <- ggplot(data@data, aes_string(x="RANK", y="P")) + 
            ggtitle(main_title) + 
            baseTheme() + 
            xlab('Rank') +
            ylab('Abundance') +
            scale_x_log10(limits=xlim,
                          breaks=scales::trans_breaks('log10', function(x) 10^x),
                          labels=scales::trans_format('log10', scales::math_format(10^.x))) +
            scale_y_continuous(labels=scales::percent) +
            geom_ribbon(aes_string(ymin="LOWER", ymax="UPPER"), fill=line_color, alpha=0.4) +
            geom_line(color=line_color)
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}


#' Plot the results of rarefyDiversity
#' 
#' \code{plotDiversityCurve} plots a \code{DiversityCurve} object.
#'
#' @param    data            \link{DiversityCurve} object returned by 
#'                           \link{rarefyDiversity}.
#' @param    colors          named character vector whose names are values in the 
#'                           \code{group} column of the \code{data} slot of \code{data},
#'                           and whose values are colors to assign to those group names.
#' @param    main_title      string specifying the plot title.
#' @param    legend_title    string specifying the legend title.
#' @param    log_x           if \code{TRUE} then plot \eqn{q} on a log scale;
#'                           if \code{FALSE} plot on a linear scale.
#' @param    log_y           if \code{TRUE} then plot the diversity/evenness scores 
#'                           on a log scale; if \code{FALSE} plot on a linear scale.
#' @param    xlim            numeric vector of two values specifying the 
#'                           \code{c(lower, upper)} x-axis limits.
#' @param    ylim            numeric vector of two values specifying the 
#'                           \code{c(lower, upper)} y-axis limits.
#' @param    annotate        string defining whether to added values to the group labels 
#'                           of the legend. When \code{"none"} (default) is specified no
#'                           annotations are added. Specifying (\code{"depth"}) adds 
#'                           sequence counts to the labels.
#' @param    score           one of \code{"diversity"} or \code{"evenness"} specifying which
#'                           score to plot on the y-asis.
#' @param    silent          if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                           object; if \code{FALSE} draw the plot.
#' @param    ...             additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  See \link{rarefyDiversity} for generating \link{DiversityCurve}
#'           objects for input. Plotting is performed with \link{ggplot}.
#' 
#' @examples
#' # Calculate diversity
#' div <- rarefyDiversity(ExampleDb, "SAMPLE", nboot=100)
#' 
#' # Plot diversity
#' plotDiversityCurve(div, legend_title="Sample")
#'
#' #' # Plot diversity
#' plotDiversityCurve(div, legend_title="Sample", score="evenness")
#' 
#' @export
plotDiversityCurve <- function(data, colors=NULL, main_title="Diversity", 
                               legend_title="Group", log_x=FALSE, log_y=FALSE,
                               xlim=NULL, ylim=NULL, annotate=c("none", "depth"), 
                               score=c("diversity", "evenness"),
                               silent=FALSE, ...) {
    # Check arguments
    annotate <- match.arg(annotate)
    score <- match.arg(score)
    
    # Define group label annotations
    if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        group_labels <- setNames(paste0(data@groups, " (N=", data@n, ")"), 
                                 data@groups)
    }
    
    # Define y-axis scores
    if (score == "diversity") {
        y_value <- "D"
        y_min <- "D_LOWER"
        y_max <- "D_UPPER"
        y_label <- expression(''^q * D)
    } else if (score == "evenness") {
        y_value <- "E"
        y_min <- "E_LOWER"
        y_max <- "E_UPPER"
        y_label <- expression(''^q * E)
    }
    
    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL
    
    # Define base plot elements
    p1 <- ggplot(data@data, aes_string(x="Q", y=y_value, group="GROUP")) + 
        ggtitle(main_title) + 
        baseTheme() + 
        xlab('q') +
        ylab(y_label) +
        geom_ribbon(aes_string(ymin=y_min, ymax=y_max, fill="GROUP"), alpha=0.4) +
        geom_line(aes_string(color="GROUP"))
    
    # Set colors and legend
    if (!is.null(colors)) {
        p1 <- p1 + scale_color_manual(name=legend_title, labels=group_labels, values=colors) +
            scale_fill_manual(name=legend_title, labels=group_labels, values=colors)
    } else {
        p1 <- p1 + scale_color_discrete(name=legend_title, labels=group_labels) +
            scale_fill_discrete(name=legend_title, labels=group_labels)
    }
    
    # Set x-axis style
    if (log_x) {
        p1 <- p1 + scale_x_continuous(trans=scales::log2_trans(), limits=xlim,
                                      breaks=scales::trans_breaks('log2', function(x) 2^x),
                                      labels=scales::trans_format('log2', scales::math_format(2^.x)))
    } else {
        p1 <- p1 + scale_x_continuous(limits=xlim)
    }
    
    # Set y-axis style
    if (log_y) {
        p1 <- p1 + scale_y_continuous(trans=scales::log2_trans(), limits=ylim,
                                      breaks=scales::trans_breaks('log2', function(x) 2^x),
                                      labels=scales::trans_format('log2', scales::math_format(2^.x)))
    } else {
        p1 <- p1 + scale_y_continuous(limits=ylim)
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))

    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}


#' Plot the results of TestDiversity
#' 
#' \code{plotDiversityTest} plots a \code{DiversityTest} object as the mean
#' with a line range indicating plus/minus one standard deviation.
#'
#' @param    data            \link{DiversityTest} object returned by 
#'                           \link{testDiversity}.
#' @param    colors          named character vector whose names are values in the 
#'                           \code{group} column of the \code{data} slot of \code{data},
#'                           and whose values are colors to assign to those group names.
#' @param    main_title      string specifying the plot title.
#' @param    legend_title    string specifying the legend title.
#' @param    log_d           if \code{TRUE} then plot the diversity scores \eqn{D} 
#'                           on a log scale; if \code{FALSE} plot on a linear scale.
#' @param    annotate        string defining whether to added values to the group labels 
#'                           of the legend. When \code{"none"} (default) is specified no
#'                           annotations are added. Specifying (\code{"depth"}) adds 
#'                           sequence counts to the labels.
#' @param    silent          if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                           object; if \code{FALSE} draw the plot.
#' @param    ...             additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  See \link{testDiversity} for generating \link{DiversityTest}
#'           objects for input. Plotting is performed with \link{ggplot}.
#' 
#' @examples
#' # All groups pass default minimum sampling threshold of 10 sequences
#' div <- testDiversity(ExampleDb, "SAMPLE", q=0, nboot=100)
#' plotDiversityTest(div, legend_title="Sample")
#' 
#' @export
plotDiversityTest <- function(data, colors=NULL, main_title="Diversity", 
                              legend_title="Group", log_d=FALSE, 
                              annotate=c("none", "depth"),
                              silent=FALSE, ...) {
    # Check arguments
    annotate <- match.arg(annotate)
    
    # Define group label annotations
    if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        group_labels <- setNames(paste0(data@groups, " (N=", data@n, ")"), 
                                 data@groups)
    }
    
    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL

    # Define plot values
    df <- data@summary %>%
        dplyr::mutate_(LOWER=~MEAN-SD, UPPER=~MEAN+SD)
    
    # Define base plot elements
    p1 <- ggplot(df, aes_string(x="GROUP")) + 
        ggtitle(main_title) + 
        baseTheme() + 
        xlab("") +
        ylab(bquote("Mean " ^ .(data@q) * D %+-% "SD")) +
        geom_linerange(aes_string(ymin="LOWER", ymax="UPPER", color="GROUP"), alpha=0.8) +
        geom_point(aes_string(y="MEAN", color="GROUP"))
    
    # Set colors and legend
    if (!is.null(colors)) {
        p1 <- p1 + scale_color_manual(name=legend_title, labels=group_labels, values=colors)
    } else {
        p1 <- p1 + scale_color_discrete(name=legend_title, labels=group_labels)
    }

    # Set x-axis style
    if (log_d) {
        p1 <- p1 + scale_y_continuous(trans=scales::log2_trans(),
                                      breaks=scales::trans_breaks('log2', function(x) 2^x),
                                      labels=scales::trans_format('log2', scales::math_format(2^.x)))
    } else {
        p1 <- p1 + scale_y_continuous()
    }

    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}