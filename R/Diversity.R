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
#' @param    group  character vector defining \code{data} columns containing grouping 
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
#' clones <- countClones(ExampleDb, group="SAMPLE")
#'
#' # With copy numbers and multiple groups
#' clones <- countClones(ExampleDb, group=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")
#' 
#' @export
countClones <- function(data, group=NULL, copy=NULL, clone="CLONE") {
    # Check input
    check <- checkColumns(data, c(clone, copy, group))
    if (check != TRUE) { stop(check) }
    
    # Tabulate clonal abundance
    if (is.null(copy)) {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize(SEQ_COUNT=n()) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(SEQ_COUNT)") %>%
            dplyr::rename_(.dots=c("CLONE"=clone))
    } else {
        clone_tab <- data %>% 
            group_by_(.dots=c(group, clone)) %>%
            dplyr::summarize_(SEQ_COUNT=interp(~length(x), x=as.name(clone)),
                              COPY_COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy))) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT")),
                           COPY_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("COPY_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(COPY_COUNT)") %>%
            dplyr::rename_(.dots=c("CLONE"=clone))
    }
    
    return(clone_tab)
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
calcInferredDiversity <- function(p, q) {
    # Use Chao estimators for unseen abundance correction
    p <- c(p, inferUnseenAbundance(p))
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

#' Creates diversity object by bootstrapping clonal distributions under rarefaction
#'
#' \code{createDiversityObject} divides a set of clones by a group annotation,
#' and resamples the sequences from each group, returning distributions of clone sizes.
#'
#' @param    data      data.frame with Change-O style columns containing clonal assignments.
#' @param    group     name of the \code{data} column containing group identifiers.
#' @param    clone     name of the \code{data} column containing clone identifiers.
#' @param    copy      name of the \code{data} column containing copy numbers for each 
#'                     sequence. If \code{copy=NULL} (the default), then clone abundance
#'                     is determined by the number of sequences. If a \code{copy} column
#'                     is specified, then clone abundances is determined by the sum of 
#'                     copy numbers within each clonal group.
#' @param    ndepth    number of observations to sample. If \code{NULL} then 
#'                     no rarefaction is performed. 
#' @param    ci        confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot     number of bootstrap realizations to generate.
#'
#' @return   A matrix of bootstrap iterations from rmultinom.

#' @export

estimateAbundance <- function(data, group, 
    clone="CLONE", copy = NULL, 
    uniform=TRUE, ndepth = NULL, nboot = 200, min_n = 30, max_n=NULL, ci=0.95){
     
	# Set confidence interval for plotting
	ci_z <- ci + (1 - ci) / 2
		
    # Private bootstrap function
    bootstrap_ <- function(data, ndepth=ndepth){
    
        # Tabulate clones
        clone_tab <- data %>% countClones(copy=copy, clone=clone)

        # Set ndepth
        ndepth_check <- clone_tab %>%
            dplyr::summarize_(SEQUENCES=interp(~sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT"))) %>%
            min(.$SEQUENCES)
        if(is.null(ndepth)){ndepth <- ndepth_check}

        # Infer Abundance
        abund_obs <- clone_tab %>% 
            select(one_of(c(clone, "SEQ_COUNT"))) %>% 
            tibble::deframe() %>%
            adjustObservedAbundance()

        # Generate bootstrap distributions
        sample_mat <- rmultinom(n=nboot, size=ndepth, abund_obs)

        return(sample_mat)
    }
    # Check input
    if (!is.data.frame(data)) {
        stop("Input data is not a data.frame")
    }
    
    # Check columns that are reported are real columns (can be NULL)
    check <- checkColumns(data, c(clone, copy, group))
    if (check != TRUE) { stop(check) }

    # Check the smallest sample depth
    sequence_ndepths <- data %>% 
        countClones(copy=copy, clone=clone, group=group) %>%
        dplyr::summarize_(SEQUENCES=interp(~sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT"))) %>%
        dplyr::select(SEQUENCES) %>% unlist() 
    
    # Check the smallest sample is not smaller than min_n
    if(min(sequence_ndepths) < min_n){stop("too few sequences in some samples")}
	
	# If max_n is set
	if(!is.null(max_n)){
		if(max(sequence_ndepths) > max_n){stop("too many sequences in some samples")}
		if(!is.null(ndepth)){if(uniform & ndepth > max_n){ndepth <- max_n}}
	}
	
    # If rarefaction is turned on, set an ndepth. Otherwise NULL ndepth (ignore).
    if(uniform & is.null(ndepth)){ndepth <- min(sequence_ndepths)} else {ndepth <- NULL}
	
    # Bootstrap abundance curves
    boot_output <- data %>%
        dplyr::group_by_(.dots=c(group)) %>%
        dplyr::do(data.frame(bootstrap_(., ndepth=ndepth)) %>% 
        tibble::rownames_to_column(clone))
    
    # Groups
    groups <- unique(data[[group]])
	
	# Create abundance slot
    abund <- boot_output %>%
        tidyr::gather(key = "N", value = "C", -one_of(c(clone, group))) %>%
        dplyr::group_by_(.dots=c(clone, group)) %>%
        dplyr::summarize(P_ERROR = qnorm(ci_z) * sd(C/ndepth), P = mean(C/ndepth)) %>%
        dplyr::mutate(LOWER = pmax(P - P_ERROR, 0), UPPER = P + P_ERROR) %>%
        dplyr::group_by_(.dots=c(group)) %>%
        dplyr::mutate(RANK = rank(-P, ties.method = "first"))

    # Create a new diversity object with bootstrap
    boot_obj <- new("AbundanceCurve",
             bootstrap=boot_output, 
			 abund=abund,
             group=group,
             groups=groups,
             clone=clone,
             uniform=uniform,
             ndepth=ndepth, 
			 nboot=nboot, 
			 ci=ci,
             min_n=min_n
    )
    
    return(boot_obj)
}


#' Generate a clonal diversity index curve
#'
#' \code{generalizeDiversity} divides a set of clones by a group annotation,
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


# Define helper functions




#' Helper function for computing alpha diversity from bootrstrap outputs
#'
#' \code{helperAlpha} divides a set of bootstrapped clones by group annotation,
#' and computes the diversity of each set. 
#'
#' @param
helperAlpha <- function(boot_output, q, clone=NULL, group=NULL){

    # Compute diversity from a column of each bootstrap
    output <- boot_output %>% 
        dplyr::ungroup() %>%
        dplyr::select(-one_of(c(clone, group))) %>%
        as.matrix() %>% 
        apply(2, calcInferredDiversity, q = q) %>%
        data.frame() %>% mutate(Q = q)

    return(output)
}


#' Helper function for computing beta diversity from bootrstrap outputs
#'
#' \code{helperBeta} divides a set of bootstrapped clones by group annotation,
#' and computes the alpha diversity. Group annotations are then ignored and 
#' gamma diversity is computed. A multiplicative beta diversity is used corresponding
#' to the gamma diversity divided by the average alpha diversity of each group.

helperBeta <- function(boot_output, q, ci_z, clone=NULL, group=NULL){

    # Compute gamma diversity metrics
    gamma <- boot_output %>%
        dplyr::group_by_(.dots=c(clone)) %>%
        dplyr::select(-one_of(c(group))) %>%
        dplyr::summarize_all(sum) %>%
        dplyr::do(helperAlpha(., q = q, clone=clone)) %>%
        tidyr::gather(key = "N", value = "GAMMA", -Q) %>%
        dplyr::mutate(GAMMA = as.numeric(GAMMA))

    # Compute alpha diversity metrics
    alpha <- boot_output %>%
        dplyr::group_by_(.dots=c(group)) %>%
        dplyr::do(helperAlpha(., q = q, clone=clone, group=group)) %>%
        dplyr::group_by(Q) %>%
        dplyr::select(-one_of(c(group))) %>%
        dplyr::summarize_all(mean) %>%
        tidyr::gather(key = "N", value = "ALPHA", -Q) %>%
        dplyr::mutate(ALPHA = as.numeric(ALPHA))

    # Perform comparisons of alpha and gamma to extract beta
    div <- bind_cols(gamma, alpha) %>%
        dplyr::group_by(Q) %>%
        dplyr::mutate(D = GAMMA/ALPHA) %>%
        dplyr::summarize(D_ERROR = qnorm(ci_z) * sd(D), D = mean(D)) %>%
        dplyr::mutate(D_LOWER = pmax(D - D_ERROR, 0), D_UPPER = D + D_ERROR)

    return(div)
}

#' Helper function for computing statistical significance
#'
#' \code{helperTest} computes the pairwise statistical significance of differences
#' in bootstrapped diversity values between two sets defined by the group column. 
#' A p-value is computed using the ECDF distribution as the frequency of bootstrap iterations
#' for which no difference is observed. 


helperTest <- function(div_df, group, q = q){

    group_pairs <- combn(unique(div_df[[group]]), 2, simplify=F)
	
    pvalue_list <- list()
    for (group_pair in group_pairs) {
        pair_list <- list()
        for(q_i in q){
            
            # Currently just testing for one diversity order
            mat1 <- div_df %>%
                dplyr::filter(.[[group]] == group_pair[1], Q == q_i) %>%
                dplyr::select(-one_of(c(group, "Q"))) %>% unlist()
            mat2 <- div_df %>%
                dplyr::filter(.[[group]] == group_pair[2], Q == q_i) %>%
                dplyr::select(-one_of(c(group, "Q"))) %>% unlist()

            if (mean(mat1) >= mean(mat2)) { g_delta <- mat1 - mat2
            } else { g_delta <- mat2 - mat1 }  

            # Compute p-value from ecdf
            p <- ecdf(g_delta)(0)
            p <- ifelse(p <= 0.5, p * 2, (1 - p) * 2)
            
            pair_list[[as.character(q_i)]] <- 
                list(DELTA_MEAN = mean(g_delta), DELTA_SD = sd(g_delta), PVALUE = p)
            
        }
        pvalue_list[[paste(group_pair, collapse=" != ")]] <- bind_rows(pair_list, .id="Q")

    }
    test_df <- bind_rows(pvalue_list, .id = "test_name")

    summary_df <- div_df %>%
        tidyr::gather(key = "N", value = "D", -one_of(c(group, "Q"))) %>%
        dplyr::mutate(D = as.numeric(D)) %>%
        dplyr::group_by_(.dots=c(group, "Q")) %>%
        dplyr::summarize(SD = sd(D), MEAN = mean(D))
	
    test_div = list(
        tests=test_df,
        summary=summary_df
    )

    return(test_div)
}


#' Function for computing the alpha diversity of bootstrap iterations
#'
#' \code{calculateAlphaDiversity} takes in a diversity object containing bootstrapped 
#' and potentially rarefied samples of clones and computes an alpha diversity metric (Hill) 
#' across a range of values (Hill values). An abundance curve is also provided along with 
#' statistical testing for significant pairwise differences in diversity. 
#' @export

calculateAlphaDiversity <- function(boot_obj, 
            min_q=0, max_q=4, step_q=0.1, 
            ci = 0.95){
    
    # Set diversity orders and confidence interval
    ci_z <- ci + (1 - ci) / 2
    q <- seq(min_q, max_q, step_q)
	if (!(0 %in% q)) { q <- c(0, q) }
		
    # Compute diversity metric for bootstrap instances
    div_df <- boot_obj@bootstrap %>%
        dplyr::group_by_(.dots=c(boot_obj@group)) %>%
        dplyr::do(helperAlpha(., q = q, clone=boot_obj@clone, group=boot_obj@group)) %>%
        dplyr::ungroup()
    
    # Summarize diversity
    div <- div_df %>%
        tidyr::gather(key = "N", value = "D", -one_of(c(boot_obj@group, "Q"))) %>%
        dplyr::mutate(D = as.numeric(D)) %>%
        dplyr::group_by_(.dots=c(boot_obj@group, "Q")) %>%
        dplyr::summarize(D_ERROR = qnorm(ci_z) * sd(D), D = mean(D)) %>%
        dplyr::mutate(D_LOWER = pmax(D - D_ERROR, 0), D_UPPER = D + D_ERROR)

    # Alpha groups
    div_groups <- unique(div[[boot_obj@group]])
    
    # Compute evenness
    div_qi <- div %>%
        filter(Q == 0) %>%
        select(one_of(c(boot_obj@group, "D")))

    div <- div %>%
        right_join(div_qi, by = boot_obj@group, suffix = c("", "_0")) %>%
        mutate(E = D/D_0, E_LOWER = D_LOWER/D_0, E_UPPER = D_UPPER/D_0) %>%
        select(-D_0)
    
    # Test
    test <- helperTest(div_df, group=boot_obj@group, q = q)
    
    boot_obj <- new("DiversityCurve",
             div=div, 
			 method="alpha",
             test=test,
             div_group=boot_obj@group,
			 div_groups=div_groups,
             q=q,  
             ci=ci)
        
   return(boot_obj) 
}

#' Function for computing the pairwise beta diversity of bootstrap iterations
#'
#' \code{calculateBetaDiversity} takes in a diversity object containing bootstrapped 
#' and potentially rarefied samples of clones and computes pairwise beta diversity metrics (Hill) 
#' across a range of values (Hill values). This is performed by using the \code{helperBeta} function
#' for pairs of groups of clones. Currently, an abundance curve is not provided as an output.  
#' Statistical testing for significant pairwise differences in diversity is also not implemented. 
#' @export

calculateBetaDiversity <- function(boot_obj, comparisons,
            min_q=0, max_q=4, step_q=0.1, 
            ci = 0.95){
    boot_obj
    # Set diversity orders and confidence interval
    ci_z <- ci + (1 - ci) / 2
    q <- seq(min_q, max_q, step_q)
	if (!(0 %in% q)) { q <- c(0, q) }
		
    # Compute pairwise beta diversity for bootstrap instances
    beta_diversity_list <- list()

    for(comparison in names(comparisons)){
        beta_diversity_list[[comparison]] <-  boot_obj@bootstrap %>%
            dplyr::ungroup() %>%
            dplyr::filter(.[[boot_obj@group]] %in% comparisons[[comparison]]) %>%
            dplyr::do(helperBeta(., q = q, clone=boot_obj@clone, group=boot_obj@group, ci_z=ci_z))
    }

    # Generate summary diversity output
    div <- bind_rows(beta_diversity_list, .id = "COMPARISON")
    
    # Beta groups
    div_groups=unique(div[["COMPARISON"]])
    
    # Compute evenness
    div_qi <- div %>%
        filter(Q == 0) %>%
        select(one_of(c("COMPARISON", "D")))

    div <- div %>%
        right_join(div_qi, by = "COMPARISON", suffix = c("", "_0")) %>%
        mutate(E = D/D_0, E_LOWER = D_LOWER/D_0, E_UPPER = D_UPPER/D_0) %>%
        select(-D_0)
		
    # Test
    #test <- testDiversity_(div_df, group=boot_obj@group, q = q)
    
    div_obj <- new("DiversityCurve",
             div=div, 
			 method="beta",
             test=NULL,
             div_group="COMPARISON",
			 div_groups=div_groups,
             q=q,  
             ci=ci)
        
	return(div_obj)
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

	bootstrap_obj <- estimateAbundance(data, group=group, clone=clone, copy=copy, nboot=nboot, min_n=min_n, max_n=max_n, uniform=uniform, ci=ci)
	diversity_obj <- calculateAlphaDiversity(bootstrap_obj, ci=ci, min_q=min_q, max_q=max_q, step_q)
	
	return(diversity_obj)	
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
                          min_n=30, max_n=NULL, nboot=2000, progress=FALSE, ci=0.95) {

	diversity_obj <- bootstrapDiversity(data, group=group, clone=clone, copy=copy, nboot=nboot, min_n=min_n, max_n=max_n, ci=ci)
	diversity_obj <- calculateAlphaDiversity(diversity_obj, min_q=q, max_q=q, step_q=0, ci=ci)
	
	return(diversity_obj)	
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


# TODO: eliminate this block


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
    
    # Check if abundance is in data
    if (is.null(data@abund)) { stop("missing abundance data") }
    
    # Check arguments
    annotate <- match.arg(annotate)
    
    # Define group label annotations
    if (all(is.na(data@groups))) {
        group_labels <- NA  
    } else if (annotate == "none") {
        group_labels <- setNames(data@groups, data@groups)
    } else if (annotate == "depth") {
        stop("option not implemented")
#         group_labels <- setNames(paste0(data$groups, " (N=", data$ndepth, ")"), 
#                                  data$groups)
    }
    
    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL
    
    if (any(!is.na(data@groups))) {
        # Define grouped plot
        p1 <- ggplot(data@abund, aes_string(x="RANK", y="P", group=data@group)) + 
            ggtitle(main_title) + 
            baseTheme() + 
            xlab('Rank') +
            ylab('Abundance') +
            scale_x_log10(limits=xlim,
                          breaks=scales::trans_breaks('log10', function(x) 10^x),
                          labels=scales::trans_format('log10', scales::math_format(10^.x))) +
            scale_y_continuous(labels=scales::percent) +
            geom_ribbon(aes_string(ymin="LOWER", ymax="UPPER", fill=data@group), alpha=0.4) +
            geom_line(aes_string(color=data@group))
        
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
        p1 <- ggplot(data@abund, aes_string(x="RANK", y="P")) + 
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
        group_labels <- setNames(data@div_groups, data@div_groups)
    } else if (annotate == "depth") {
        stop("option not implemented")
#         group_labels <- setNames(paste0(data$groups, " (N=", data$ndepth, ")"), 
#                                  data$groups)
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
    p1 <- ggplot(data@div, aes_string(x="Q", y=y_value, group=data@div_group)) + 
        ggtitle(main_title) + 
        baseTheme() + 
        xlab('q') +
        ylab(y_label) +
        geom_ribbon(aes_string(ymin=y_min, ymax=y_max, fill=data@div_group), alpha=0.4) +
        geom_line(aes_string(color=data@div_group))
    
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
#' @param    q_i			 diversity index at which to display differences. 
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
plotDiversityTest <- function(data, q_i=NULL, colors=NULL, main_title="Diversity", 
                              legend_title="Group", log_d=FALSE, 
                              annotate=c("none", "depth"),
                              silent=FALSE, ...) {
    
    # Check if abundance is in data
    if (is.null(data@test)) { stop("missing testing from data") }
    
	# Set q_i if NULL (when using testDiversity function)
	if(is.null(q_i)){q_i <- unique(unlist(data@test$test["Q"]))[-1]}
	
    # Check if q is in data
    if (!(q_i %in% data@q)) { stop("hill index not assessed") }
    
    # Check arguments
    annotate <- match.arg(annotate)
    
    # Define group label annotations
    if (annotate == "none") {
        group_labels <- setNames(data@div_groups, data@div_groups)
    } else if (annotate == "depth") {
        stop("option not implemented")
#         group_labels <- setNames(paste0(data$groups, " (N=", data$ndepth, ")"), 
#                                  data$groups)
    }
    
    # Stupid hack for check NOTE about `.x` in math_format
    .x <- NULL

    # Define plot values
    df <- data@test$summary %>%
        dplyr::filter(Q == q_i) %>%
        dplyr::mutate_(LOWER=~MEAN-SD, UPPER=~MEAN+SD)
    
    # Define base plot elements
    p1 <- ggplot(df, aes_string(x=data@div_group)) + 
        ggtitle(main_title) + 
        baseTheme() + 
        xlab("") +
        ylab(bquote("Mean " ^ .(data@q) * D %+-% "SD")) +
        geom_linerange(aes_string(ymin="LOWER", ymax="UPPER", color=data@div_group), alpha=0.8) +
        geom_point(aes_string(y="MEAN", color=data@div_group))
    
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