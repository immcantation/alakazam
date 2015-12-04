# Clonality analysis

#' @include Alakazam.R
NULL

#### Classes ####


#### Methods ####


#### Calculation functions ####

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
#' @return  A data.frame summarizing clone counts and frequencies with columns:
#'          \itemize{
#'            \item \code{clone}:       clone identifier.
#'            \item \code{seq_count}:   total number of sequences for the clone.
#'            \item \code{seq_freq}:    frequency of the clone as a fraction of the total
#'                                      number of sequences within each group.
#'            \item \code{copy_count}:  sum of the copy counts in the \code{copy} column.
#'                                      Only present if the \code{copy} argument is 
#'                                      specified.
#'            \item \code{copy_freq}:   frequency of the clone as a fraction of the total
#'                                      copy number within each group. Only present if 
#'                                      the \code{copy} argument is specified.
#'          }
#'          Also includes sdditional columns specified in the \code{groups} argument.
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Without copy numbers
#' clones <- countClones(df, groups="SAMPLE")
#'
#' # With copy numbers and multiple groups
#' clones <- countClones(df, groups=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")
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
            dplyr::summarize(seq_count=n()) %>%
            dplyr::mutate(seq_freq=seq_count/sum(seq_count, na.rm=TRUE)) %>%
            dplyr::arrange(desc(seq_count)) %>%
            dplyr::rename_(.dots=c("clone"=clone))
    } else {
        clone_tab <- data %>% 
            group_by_(.dots=c(groups, clone)) %>%
            dplyr::summarize_(seq_count=interp(~length(x), x=as.name(clone)),
                              copy_count=interp(~sum(x, na.rm=TRUE), x=as.name(copy))) %>%
            dplyr::mutate(seq_freq=seq_count/sum(seq_count, na.rm=TRUE),
                          copy_freq=copy_count/sum(copy_count, na.rm=TRUE)) %>%
            dplyr::arrange(desc(copy_count)) %>%
            dplyr::rename_(.dots=c("clone"=clone))
    }

    return(clone_tab)
}


#' Estimates the complete clonal relative abundance distribution
#' 
#' \code{estimateAbundance} estimates the complete clonal relative abundance distribution 
#' and confidence intervals on clone sizes using bootstrapping.
#' 
#' @param    data   data.frame with Change-O style columns containing clonal assignments.
#' @param    group  name of the \code{data} column containing group identifiers.
#' @param    clone  name of the \code{data} column containing clone identifiers.
#' @param    copy   name of the \code{data} column containing copy numbers for each 
#'                  sequence. If \code{copy=NULL} (the default), then clone abundance
#'                  is determined by the number of sequences. If a \code{copy} column
#'                  is specified, then clone abundances is determined by the sum of 
#'                  copy numbers within each clonal group.
#' @param    ci     confidence interval to calculate; the value must be between 0 and 1.
#' @param    nboot  number of bootstrap realizations to generate.
#' 
#' @return   A data.frame with relative clonal abundance data and confidence intervals,
#'           containing the following columns:
#'           \itemize{
#'             \item  \code{group}:  group identifier.
#'             \item  \code{clone}:  clone identifier.
#'             \item  \code{p}:      relative abundance of the clone.
#'             \item  \code{lower}:  lower confidence inverval bound.
#'             \item  \code{upper}:  upper confidence interval bound.
#'             \item  \code{rank}:   the rank of the clone abundance.
#'           }
#'           
#' @details
#' The complete clonal abundance distribution determined inferred by using the Chao1 
#' estimator to estimate the number of seen clones, and applying the relative abundance 
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
#' See \code{\link{rarefyDiversity}} for a similar application to clonal diversity.
#'           
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' abund <- estimateAbundance(df, "SAMPLE", nboot=100)
#'
#' @export
estimateAbundance <- function(data, group, clone="CLONE", copy=NULL, ci=0.95, nboot=2000) {
    #group="SAMPLE"; clone="CLONE"; copy="UID_CLUSTCOUNT"; ci=0.95; nboot=200

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
    
    # Summarize groups
    group_tab <- clone_tab %>%
        group_by_(.dots=c(group)) %>%
        dplyr::summarize(SEQUENCES=sum(COUNT, na.rm=TRUE))
    group_set <- as.character(group_tab[[group]])
    
    # Set confidence interval
    ci_z <- ci + (1 - ci) / 2
    
    # Generate diversity index and confidence intervals via resampling
    cat("-> ESTIMATING ABUNDANCE\n")
    pb <- txtProgressBar(min=0, max=length(group_set), initial=0, width=40, style=3)
    abund_list <- list()
    i <- 0
    for (g in group_set) {
        i <- i + 1
        nsam <- group_tab$SEQUENCES[group_tab[[group]] == g]
        
        # Infer complete abundance distribution
        # TODO:  can be a single function (wrapper) for both this and rarefyDiversity
        abund_obs <- clone_tab$COUNT[clone_tab[[group]] == g]
        p1 <- adjustObservedAbundance(abund_obs)
        p2 <- inferUnseenAbundance(abund_obs)
        p <- c(p1, p2)
        names(p) <- c(clone_tab$CLONE[clone_tab[[group]] == g], 
                      paste0("U", 1:length(p2)))
        
        # Bootstrap abundance
        boot_mat <- rmultinom(nboot, nsam, p) / nsam
        
        # Assign confidence intervals based on variance of bootstrap realizations
        boot_sd <- apply(boot_mat, 1, sd)
        boot_err <- qnorm(ci_z) * boot_sd
        p_lower <- pmax(p - boot_err, 0)
        p_upper <- p + boot_err
        
        # Assemble and sort abundance data.frame
        abund_df <- dplyr::data_frame(clone=names(p), p=p, lower=p_lower, upper=p_upper)
        abund_df <- dplyr::arrange(abund_df, desc(p))
        abund_df$rank <- 1:nrow(abund_df)
        abund_list[[g]] <- abund_df
        
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    
    return(bind_rows(abund_list, .id="group"))
}


#### Plotting functions ####

#' Plots a clonal abundance distribution
#' 
#' \code{plotAbundance} plots the results from estimating the complete clonal relative 
#' abundance distribution. The distribution is plotted as a log rank abundance 
#' distribution.
#' 
#' @param    data          data.frame returned by \link{estimateAbundance}.
#' @param    colors        named character vector whose names are values in the 
#'                         \code{group} column of \code{data} and whose values are 
#'                         colors to assign to those group names.
#' @param    main_title    string specifying the plot title.
#' @param    legend_title  string specifying the legend title.
#' @param    xlim          numeric vector of two values specifying the 
#'                         \code{c(lower, upper)} x-axis limits.
#' @param    ylim          numeric vector of two values specifying the 
#'                         \code{c(lower, upper)} y-axis limits.
#' @param    silent        if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                         object; if \code{FALSE} draw the plot.
#' @param    ...           additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  
#' See \code{\link{estimateAbundance}} for generating the input abundance distribution.
#' Plotting is performed with \code{\link{ggplot}}.
#'           
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Plot
#' abund <- estimateAbundance(df, "SAMPLE", nboot=100)
#' plotAbundance(abund)
#' 
#' @export
plotAbundance <- function(data, colors=NULL, main_title="Rank Abundance", 
                          legend_title=NULL, xlim=NULL, ylim=NULL, 
                          silent=FALSE, ...) {
    # TODO: additional styles. rank abundance, box/violin
    
    # Define base plot elements
    g1 <- ggplot(data, aes(x=rank, y=p, group=group)) + 
        ggtitle(main_title) + 
        getBaseTheme() + 
        xlab('Rank') +
        ylab('Abundance') +
        scale_x_log10(limits=xlim,
                      breaks=trans_breaks('log10', function(x) 10^x),
                      labels=trans_format('log10', math_format(10^.x))) +
        scale_y_continuous(labels=percent) +
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=group), alpha=0.4) +
        geom_line(aes(color=group))
    
    # Set colors and legend
    if (!is.null(colors)) {
        g1 <- g1 + scale_color_manual(name=legend_title, values=colors) +
            scale_fill_manual(name=legend_title, values=colors)
    } else {
        g1 <- g1 + scale_color_discrete(name=legend_title) +
            scale_fill_discrete(name=legend_title)
    }

    # Add additional theme elements
    g1 <- g1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(g1) }
    
    invisible(g1)
}
