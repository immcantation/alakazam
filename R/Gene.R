# Gene usage analysis

#### Calculation functions ####

#' Tabulates V(D)J allele, gene or family usage.
#' 
#' Determines the count and relative abundance of V(D)J alleles, genes or families within
#' groups.
#'
#' @param    data    data.frame with Change-O style columns.
#' @param    gene    column containing allele assignments. Only the first allele in the
#'                   column will be considered.
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
#' @param    mode    one of \code{c("gene", "family", "allele")} defining
#'                   the degree of specificity regarding allele calls. Determines whether 
#'                   to return counts for genes, families or alleles.
#' 
#' @return   A data.frame summarizing family, gene or allele counts and frequencies 
#'           with columns:
#'           \itemize{
#'             \item \code{GENE}:        name of the family, gene or allele
#'             \item \code{SEQ_COUNT}:   total number of sequences, or clones, for the gene.
#'             \item \code{SEQ_FREQ}:    frequency of the gene as a fraction of the total
#'                                       number of sequences, or clones, within each grouping.
#'             \item \code{COPY_COUNT}:  sum of the copy counts in the \code{copy} column.
#'                                       for each gene. Only present if the \code{copy} 
#'                                       argument is specified.
#'             \item \code{COPY_FREQ}:   frequency of the gene as a fraction of the total
#'                                       copy number within each group. Only present if 
#'                                       the \code{copy} argument is specified.
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
#'@export
countGenes <- function(data, gene, groups=NULL, copy=NULL, clone=NULL,
                       mode=c("gene", "allele", "family")) {
    ## DEBUG
    # data=db; gene="V_CALL"; groups=NULL; mode="gene"; clone="CLONE"
    # data=subset(db, CLONE == 3138)
    
    # Check input
    mode <- match.arg(mode)
    check <- checkColumns(data, c(gene, groups, copy))
    if (check != TRUE) { stop(check) }
    
    # Subset to one sequence per clone if required
    if (!is.null(clone) & is.null(copy)) {
        # Find count of genes within each clone
        data <- data %>%
            group_by_(.dots=c(groups, clone, gene)) %>%
            dplyr::mutate(CLONE_GENE_COUNT=n())
        # Keep first row that corresponds to the most common gene
        data <- data %>%
            group_by_(.dots=c(groups, clone)) %>%
            slice_(interp(~which.max(x), x=as.name("CLONE_GENE_COUNT"))) %>%
            ungroup() %>%
            select_(interp(~-x, x=as.name("CLONE_GENE_COUNT")))
    } else if (!is.null(clone) & !is.null(copy)) {
        warning("Specifying both 'copy' and 'clone' columns is not meaningful. ",
                "The 'clone' argument will be ignored.")
    }
    # Extract gene, allele or family assignments
    gene_func <- switch(mode,
                        allele=getAllele,
                        gene=getGene,
                        family=getFamily)
    data[[gene]] <- gene_func(data[[gene]], first=TRUE)
    
    # Tabulate clonal abundance
    if (is.null(copy)) {
        gene_tab <- data %>% 
            group_by_(.dots=c(groups, gene)) %>%
            dplyr::summarize(SEQ_COUNT=n()) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(SEQ_COUNT)") %>%
            dplyr::rename_(.dots=c("GENE"=gene))
    } else {
        gene_tab <- data %>% 
            group_by_(.dots=c(groups, gene)) %>%
            dplyr::summarize_(SEQ_COUNT=interp(~length(x), x=as.name(gene)),
                              COPY_COUNT=interp(~sum(x, na.rm=TRUE), x=as.name(copy))) %>%
            dplyr::mutate_(SEQ_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("SEQ_COUNT")),
                           COPY_FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("COPY_COUNT"))) %>%
            dplyr::arrange_(.dots="desc(COPY_COUNT)") %>%
            dplyr::rename_(.dots=c("GENE"=gene))
    }
    
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
#' @export
getSegment <- function(segment_call, segment_regex, first=TRUE, collapse=TRUE, 
                       strip_d=TRUE, sep=",") {
    # Define boundaries of individual segment calls
    edge_regex <- paste0("[^", sep, "]*")
    
    # Extract calls
    r <- gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), "\\1", 
              segment_call, perl=T)
    
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
getAllele <- function(segment_call, first=TRUE, collapse=TRUE, strip_d=TRUE, sep=",") {    
    allele_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+[-/\\w]*[-\\*][\\.\\w]+)'
    r <- getSegment(segment_call, allele_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, sep=sep)
    
    return(r)
}


#' @rdname getSegment
#' @export
getGene <- function(segment_call, first=TRUE, collapse=TRUE, strip_d=TRUE, sep=",") {
    gene_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+[-/\\w]*)'
    r <- getSegment(segment_call, gene_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, sep=sep)
    
    return(r)
}


#' @rdname getSegment
#' @export
getFamily <- function(segment_call, first=TRUE, collapse=TRUE, strip_d=TRUE, sep=",") {
    family_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+)'
    r <- getSegment(segment_call, family_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, sep=sep)
    
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
