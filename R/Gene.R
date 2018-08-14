# Gene usage analysis

#### Calculation functions ####

#' Tabulates V(D)J allele, gene or family usage.
#' 
#' Determines the count and relative abundance of V(D)J alleles, genes or families within
#' groups.
#'
#' @param    data    data.frame with Change-O style columns.
#' @param    gene    column containing allele assignments. Only the first allele in the
#'                   column will be considered when \code{mode} is "gene", "family" or 
#'                   "allele". The value will be used as it is with \code{mode="asis"}. 
#' @param    groups  columns containing grouping variables. If \code{NULL} do not group.
#' @param    copy    name of the \code{data} column containing copy numbers for each 
#'                   sequence. If this value is specified, then total copy abundance
#'                   is determined by the sum of copy numbers within each gene.
#'                   This argument is ignored if \code{clone} is specified.
#' @param    clone   name of the \code{data} column containing clone identifiers for each 
#'                   sequence. If this value is specified, then genes will be counted only
#'                   once for each clone. Note, this is accomplished by using the most 
#'                   common gene within each \code{clone} identifier. As such,
#'                   ambiguous alleles within a clone will not be accurately represented.
#' @param    mode    one of \code{c("gene", "family", "allele", "asis")} defining
#'                   the degree of specificity regarding allele calls. Determines whether 
#'                   to return counts for genes (calling \code{getGene}), 
#'                   families (calling \code{getFamily}), alleles (calling 
#'                   \code{getAllele}) or using the value as it is in the column
#'                   \code{gene}, without any processing.
#' @param    fill  logical of \code{c(TRUE, FALSE)} specifying when if groups (when specified)
#'                   lacking a particular gene should be counted as 0 if TRUE or not (omitted) 
#' 
#' @return   A data.frame summarizing family, gene or allele counts and frequencies 
#'           with columns:
#'           \itemize{
#'             \item \code{GENE}:         name of the family, gene or allele
#'             \item \code{SEQ_COUNT}:    total number of sequences for the gene.
#'             \item \code{SEQ_FREQ}:     frequency of the gene as a fraction of the total
#'                                        number of sequences within each grouping.
#'             \item \code{COPY_COUNT}:   sum of the copy counts in the \code{copy} column.
#'                                        for each gene. Only present if the \code{copy} 
#'                                        argument is specified.
#'             \item \code{COPY_FREQ}:    frequency of the gene as a fraction of the total
#'                                        copy number within each group. Only present if 
#'                                        the \code{copy} argument is specified.
#'             \item \code{CLONE_COUNT}:  total number of clones for the gene.
#'             \item \code{CLONE_FREQ}:   frequency of the gene as a fraction of the total
#'                                        number of clones within each grouping.
#'           }
#'           Additional columns defined by the \code{groups} argument will also be present.
#'
#' @examples
#' # Without copy numbers
#' genes <- countGenes(ExampleDb, gene="V_CALL", groups="SAMPLE", mode="family")
#' genes <- countGenes(ExampleDb, gene="V_CALL", groups="SAMPLE", mode="gene")
#' genes <- countGenes(ExampleDb, gene="V_CALL", groups="SAMPLE", mode="allele")
#'
#' # With copy numbers and multiple groups
#' genes <- countGenes(ExampleDb, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
#'                     copy="DUPCOUNT", mode="family")
#' 
#' # Count by clone
#' genes <- countGenes(ExampleDb, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
#'                     clone="CLONE", mode="family")
#'
#' # Count absent genes 
#' genes <- countGenes(ExampleDb, gene="V_CALL", groups="SAMPLE", mode="allele", fill = TRUE)
#'
#'@export
countGenes <- function(data, gene, groups=NULL, copy=NULL, clone=NULL, fill=FALSE,
                       mode=c("gene", "allele", "family", "asis")) {
    ## DEBUG
    # data=ExampleDb; gene="V_CALL"; groups=NULL; mode="gene"; clone="CLONE"
    # data=subset(db, CLONE == 3138)
    
    # Check input
    mode <- match.arg(mode)
    check <- checkColumns(data, c(gene, groups, copy))
    if (check != TRUE) { stop(check) }

    # Extract gene, allele or family assignments
    if (mode != "asis") {
        gene_func <- switch(mode,
                            allele=getAllele,
                            gene=getGene,
                            family=getFamily)
        data[[gene]] <- gene_func(data[[gene]], first=TRUE)
    }
    
    # Tabulate abundance
    if (is.null(copy) & is.null(clone)) {
        # Tabulate sequence abundance
        gene_tab <- data %>% 
            group_by_(.dots=c(groups, gene)) %>%
            dplyr::summarize(SEQ_COUNT=n()) %>%
            mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT"))) %>%
            arrange_(.dots="desc(SEQ_COUNT)")
    } else if (!is.null(clone) & is.null(copy)) {
        # Find count of genes within each clone and keep first with maximum count
        gene_tab <- data %>%
            group_by_(.dots=c(groups, clone, gene)) %>%
            dplyr::mutate(CLONE_GENE_COUNT=n()) %>%
            ungroup() %>%
            group_by_(.dots=c(groups, clone)) %>%
            slice_(interp(~which.max(x), x=as.name("CLONE_GENE_COUNT"))) %>%
            ungroup() %>%
            group_by_(.dots=c(groups, gene)) %>%
            dplyr::summarize(CLONE_COUNT=n()) %>%
            mutate_(CLONE_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("CLONE_COUNT"))) %>%
            arrange_(.dots="desc(CLONE_COUNT)")
    } else {
        if (!is.null(clone) & !is.null(copy)) {
            warning("Specifying both 'copy' and 'clone' columns is not meaningful. ",
                    "The 'clone' argument will be ignored.")
        }
        # Tabulate copy abundance
        gene_tab <- data %>% 
            group_by_(.dots=c(groups, gene)) %>%
            summarize_(SEQ_COUNT=interp(~length(x), x=as.name(gene)),
                       COPY_COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy))) %>%
            mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT")),
                    COPY_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("COPY_COUNT"))) %>%
            arrange_(.dots="desc(COPY_COUNT)")
    }

    # If a gene is present in one GROUP but not another, will fill the COUNT and FREQ with 0s
    if(fill){
        gene_tab <- gene_tab %>%
            ungroup() %>%
            complete_(as.list(c(groups, gene))) %>%
            replace(is.na(.), 0)
    }

    # Rename gene column
    gene_tab <- rename_(gene_tab, .dots=c("GENE"=gene))
    
    return(gene_tab)
}


#### Annotation functions ####

#' Get Ig segment allele, gene and family names
#' 
#' \code{getSegment} performs generic matching of delimited segment calls with a custom 
#' regular expression. \link{getAllele}, \link{getGene} and \link{getFamily} extract 
#' the allele, gene and family names, respectively, from a character vector of 
#' immunoglobulin (Ig) or TCR segment allele calls in IMGT format.
#'
#' @param     segment_call    character vector containing segment calls delimited by commas.
#' @param     segment_regex   string defining the segment match regular expression.
#' @param     first           if \code{TRUE} return only the first call in 
#'                            \code{segment_call}; if \code{FALSE} return all calls 
#'                            delimited by commas.
#' @param     collapse        if \code{TRUE} check for duplicates and return only unique 
#'                            segment assignments; if \code{FALSE} return all assignments 
#'                            (faster). Has no effect if \code{first=TRUE}.
#' @param     strip_d         if \code{TRUE} remove the "D" from the end of gene annotations 
#'                            (denoting a duplicate gene in the locus); 
#'                            if \code{FALSE} do not alter gene names.
#' @param     omit_nl         if \code{TRUE} remove non-localized (NL) genes from the result.
#'                            Only applies at the gene or allele level.
#' @param     sep             character defining both the input and output segment call 
#'                            delimiter.
#' 
#' @return    A character vector containing allele, gene or family names.
#' 
#' @references
#'   \url{http://imgt.org}
#'
#' @seealso  \link{countGenes}
#'
#' @examples
#' # Light chain examples
#' kappa_call <- c("Homsap IGKV1D-39*01 F,Homsap IGKV1-39*02 F,Homsap IGKV1-39*01",
#'                 "Homsap IGKJ5*01 F")
#'
#' getAllele(kappa_call)
#' getAllele(kappa_call, first=FALSE)
#' getAllele(kappa_call, first=FALSE, strip_d=FALSE)
#' 
#' getGene(kappa_call)
#' getGene(kappa_call, first=FALSE)
#' getGene(kappa_call, first=FALSE, strip_d=FALSE)
#' 
#' getFamily(kappa_call)
#' getFamily(kappa_call, first=FALSE)
#' getFamily(kappa_call, first=FALSE, collapse=FALSE)
#' getFamily(kappa_call, first=FALSE, strip_d=FALSE)
#' 
#' # Heavy chain examples
#' heavy_call <- c("Homsap IGHV1-69*01 F,Homsap IGHV1-69D*01 F", 
#'                 "Homsap IGHD1-1*01 F", 
#'                 "Homsap IGHJ1*01 F")
#' 
#' getAllele(heavy_call, first=FALSE)
#' getAllele(heavy_call, first=FALSE, strip_d=FALSE)
#'
#' getGene(heavy_call, first=FALSE)
#' getGene(heavy_call, first=FALSE, strip_d=FALSE)
#'
#' # Filtering non-localized genes
#' nl_call <- c("IGHV3-NL1*01,IGHV3-30-3*01,IGHV3-30*01", 
#'              "Homosap IGHV3-30*01 F,Homsap IGHV3-NL1*01 F",
#'              "IGHV1-NL1*01")
#'              
#' getAllele(nl_call, first=FALSE, omit_nl=TRUE)
#' getGene(nl_call, first=FALSE, omit_nl=TRUE)
#' getFamily(nl_call, first=FALSE, omit_nl=TRUE)
#'
#' @export
getSegment <- function(segment_call, segment_regex, first=TRUE, collapse=TRUE, 
                       strip_d=TRUE, omit_nl=FALSE, sep=",") {
    # Define boundaries of individual segment calls
    edge_regex <- paste0("[^", sep, "]*")
    
    # Extract calls
    r <- gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), "\\1", 
              segment_call, perl=T)
    
    # Remove NL genes
    if (omit_nl) {
        nl_regex <- paste0('(IG[HLK]|TR[ABGD])[VDJ][0-9]+-NL[0-9]([-/\\w]*[-\\*][\\.\\w]+)*(', 
                           sep, "|$)")
        r <- gsub(nl_regex, "", r, perl=TRUE)
    }
    
    # Strip D from gene names if required
    if (strip_d) {
        strip_regex <- paste0("(?<=[A-Z0-9])D(?=\\*|-|", sep, "|$)")
        r <- gsub(strip_regex, "", r, perl=TRUE)
    }
    
    # Collapse to unique set if required
    if (first) {
        r <- gsub(paste0(sep, ".*$"), "", r)
    } else if (collapse) {
        r <- sapply(strsplit(r, sep), function(x) paste(unique(x), collapse=sep))
    }
    
    return(r)
}


#' @rdname getSegment
#' @export
getAllele <- function(segment_call, first=TRUE, collapse=TRUE, 
                      strip_d=TRUE, omit_nl=FALSE, sep=",") {    
    allele_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+[-/\\w]*[-\\*][\\.\\w]+)'
    r <- getSegment(segment_call, allele_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, omit_nl=omit_nl, sep=sep)
    
    return(r)
}


#' @rdname getSegment
#' @export
getGene <- function(segment_call, first=TRUE, collapse=TRUE, 
                    strip_d=TRUE, omit_nl=FALSE, sep=",") {
    gene_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+[-/\\w]*)'
    r <- getSegment(segment_call, gene_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, omit_nl=omit_nl, sep=sep)
    
    return(r)
}


#' @rdname getSegment
#' @export
getFamily <- function(segment_call, first=TRUE, collapse=TRUE, 
                      strip_d=TRUE, omit_nl=FALSE, sep=",") {
    family_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+)'
    r <- getSegment(segment_call, family_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, omit_nl=omit_nl, sep=sep)
    
    return(r)
}


#### Utility functions ####

#' Sort V(D)J genes
#'
#' \code{sortGenes} sorts a vector of V(D)J gene names by either lexicographic ordering 
#' or locus position. 
#' 
#' @param    genes    vector of strings respresenting V(D)J gene names.
#' @param    method   string defining the method to use for sorting genes. One of:
#'                    \itemize{
#'                      \item \code{"name"}:      sort in lexicographic order. Order is by 
#'                                                family first, then gene, and then allele. 
#'                      \item \code{"position"}:  sort by position in the locus, as
#'                                                determined by the final two numbers 
#'                                                in the gene name. Non-localized genes 
#'                                                are assigned to the highest positions.
#'                    }
#'                    
#' @return   A sorted character vector of gene names.
#' 
#' @seealso  See \code{getAllele}, \code{getGene} and \code{getFamily} for parsing
#'           gene names.
#' 
#' @examples
#' # Create a list of allele names
#' genes <- c("IGHV1-69D*01","IGHV1-69*01","IGHV4-38-2*01","IGHV1-69-2*01",
#'            "IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
#'            "IGHV1-2", "IGHV1-2*02", "IGHV1-69*02")
#' 
#' # Sort genes by name
#' sortGenes(genes)
#' 
#' # Sort genes by position in the locus
#' sortGenes(genes, method="pos")
#' 
#' @export
sortGenes <- function(genes, method=c("name", "position")) { 
    ## DEBUG
    # method="name"
    
    # Check arguments
    method <- match.arg(method)

    # Build sorting table
    sort_tab <- data_frame(CALL=sort(getAllele(genes, first=FALSE, strip_d=FALSE))) %>%
        # Determine the gene and family
        mutate_(FAMILY=interp(~getFamily(x, first=TRUE, strip_d=FALSE), x=as.name("CALL")),
                GENE=interp(~getGene(x, first=TRUE, strip_d=FALSE), x=as.name("CALL")),
                ALLELE=interp(~getAllele(x, first=TRUE, strip_d=FALSE), x=as.name("CALL"))) %>%
        # Identify first gene number, second gene number and allele number
        mutate_(G1=interp(~gsub("[^-]+-([^-\\*D]+).*", "\\1", x), x=as.name("GENE")),
                G1=interp(~as.numeric(gsub("[^0-9]+", "99", x)), x=as.name("G1")),
                G2=interp(~gsub("[^-]+-[^-]+-?", "", x), x=as.name("GENE")),
                G2=interp(~as.numeric(gsub("[^0-9]+", "99", x)), x=as.name("G2")),
                A1=interp(~as.numeric(sub("[^\\*]+\\*|[^\\*]+$", "", x)), x=as.name("ALLELE")))

    # Convert missing values to 0
    sort_tab[is.na(sort_tab)] <- 0
    
    # Sort
    if (method == "name") {  
        sorted_genes <- arrange_(sort_tab, ~FAMILY, ~G1, ~G2, ~A1)[["CALL"]]
    } else if (method == "position") {
        sorted_genes <- arrange_(sort_tab, ~desc(G1), ~desc(G2), ~FAMILY, ~A1)[["CALL"]]
    }
    
    return(sorted_genes)
}
