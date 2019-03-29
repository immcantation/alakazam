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
#' genes <- countGenes(ExampleDb, gene="V_CALL", groups="SAMPLE", 
#'                     mode="allele", fill=TRUE)
#'
#'@export
countGenes <- function(data, gene, groups=NULL, copy=NULL, clone=NULL, fill=FALSE,
                       mode=c("gene", "allele", "family", "asis")) {
    ## DEBUG
    # data=ExampleDb; gene="V_CALL"; groups=NULL; mode="gene"; clone="CLONE"
    # data=subset(db, CLONE == 3138)
    # Hack for visibility of dplyr variables
    . <- NULL
    
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
    if (fill) {
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
    allele_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+[-/\\w]*[-\\*]*[\\.\\w]+)'
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




#' Group sequences by gene assignment
#'
#' \code{groupGenes} will group rows by shared V and J gene assignments, 
#' and optionally also by junction lengths.
#' Both VH:VL paired single-cell BCR-seq and unpaired bulk-seq (heavy chain-only)
#' are supported.
#' In the case of ambiguous (multiple) gene assignments, the grouping may
#' be specified to be a union across all ambiguous V and J gene pairs, 
#' analagous to single-linkage clustering (i.e., allowing for chaining).
#'
#' @param    data    data.frame containing sequence data.
#' @param    v_call  name of the column containing the heavy chain V-segment allele calls.
#' @param    j_call  name of the column containing the heavy chain J-segment allele calls.
#' @param    first   if \code{TRUE} only the first call of the gene assignments 
#'                   is used. if \code{FALSE} the union of ambiguous gene 
#'                   assignments is used to group all sequences with any 
#'                   overlapping gene calls.
#' @param    separator_within_seq    a single character specifying the separator
#'                                   between multiple annotations for a single
#'                                   sequence. Defaults to \code{,}.
#' @param    separator_between_seq   a single character specifying the separator
#'                                   between multiple sequences. Defaults to \code{;}.
#' @param    single_cell_with_light  a Boolean value specifying whether single-cell
#'                                   mode with VH:VL paired BCR input is being supplied.
#'                                   Defaults to \code{FALSE}.
#' @param    v_call_light            name of the column containing the light chain V-segment
#'                                   allele calls. Only applicable for single-cell mode.
#' @param    j_call_light            name of the column containing the light chain J-segment
#'                                   allele calls. Only applicable for single-cell mode.
#' @param    junc_len_heavy          name of the column containing the heavy chain junction
#'                                   length. Optional.
#' @param    junc_len_light          name of the column containing the light chain junction
#'                                   length. Optional and only applicable for single-cell mode.
#'
#' @return   Returns a modified \code{data} data.frame with union indices 
#'           in the \code{VJ_GROUP} column.
#'
#' @details
#' All rows containing \code{NA} valies in their \code{v_call}, \code{j_call},
#' \code{v_call_light}, \code{j_call_light}, \code{junc_len_heavy}, \code{junc_len_light} 
#' (if specified) columns will be removed. 
#' A warning will be issued when a row containing an \code{NA} is removed.
#' 
#' Ambiguous gene assignments are assumed to be separated by commas.
#' 
#' @section Expectation for single-cell input with VH:VL pairing
#' With \code{single_cell_with_light=TRUE}, it is assumed that 
#'   \itemize{
#'      \item every row represents a single cell
#'      \item each cell possibly contains multiple heavy and/pr light chains
#'      \item multiple chains, if any, are separated by \code{separator_between_seq}
#'            without any space
#'   }
#'   
#' Every chain must have its own V(D)J annotation. Example:
#' 
#' A cell/row has 1 heavy chain and 2 light chains.
#' 
#' The two light chains are separated by \code{;}.   
#' 
#' V annotations for the light chains: \code{Homsap IGKV1-39*01 F,Homsap IGKV1D-39*01 F;Homsap IGLV2-11*01 F}.
#' 
#' J annotations for the light chains: \code{Homsap IGLJ2*01 F;Homsap IGLJ3*01 F}.
#' 
#' Notice that for both the V and the J annotations, there is a \code{;} separating the 
#' two annotations for the two light chains.
#' 
#' It cannot be the case that there are 2 V annotations but only 1 J annotation for 
#' the two light chains. Both J annotations must be spelled out for each light chain,
#' separated by \code{separator_between_seq}, even if the annotated alleles are the same.
#' 
#' This one-to-one annotation-to-chain correspondence for both V and J is explicitly
#' checked and an error raised if the requirement is not met. 
#' 
#' @examples
#' # Group by genes
#' db <- groupGenes(data=ExampleDb)
#'  
#' @export
groupGenes <- function(data, v_call="V_CALL", j_call="J_CALL", first=FALSE,
                       separator_within_seq=",", separator_between_seq=";", 
                       single_cell_with_light=F, v_call_light, j_call_light, 
                       junc_len_heavy=NULL, junc_len_light=NULL) {
    
    # if using data in which each row represents a cell with both heavy and light chains
    if (single_cell_with_light) {
        
        # argument check for `v_call_light` and `j_call_light`
        # must be specified
        missing(v_call_light)
        missing(j_call_light)
        
        # must be in colnames of `data`
        if (!all( c(v_call_light, j_call_light) %in% colnames(data) )) {
            stop("v_call_light and/or j_call_light not found as a column in data")
        }
        
        # one-to-one annotation-to-chain correspondence for both V and J (light)
        # for each cell/row, number of bewteen_seq separators in light V annotation and in light J annotation must match
        n_separator_btw_seq_v_light <- stringi::stri_count_fixed(str=data[[v_call_light]], pattern=separator_between_seq)
        n_separator_btw_seq_j_light <- stringi::stri_count_fixed(str=data[[j_call_light]], pattern=separator_between_seq)
        if (!all( n_separator_btw_seq_v_light == n_separator_btw_seq_j_light )) {
            stop("Requirement not met: one-to-one annotation-to-chain correspondence for both V and J (light)")
        }
        
    } else {
        v_call_light <- NULL
        j_call_light <- NULL
    }
    
    # one-to-one annotation-to-chain correspondence for both V and J (heavy)
    # for each cell/row, number of bewteen_seq separators in heavy V annotation and in heavy J annotation must match
    # (in theory, there should be 1 heavy chain per cell; but 10x can return cell with >1 heavy chains and 
    #  you never know if the user will supply this cell as input)
    n_separator_btw_seq_v_heavy <- stringi::stri_count_fixed(str=data[[v_call]], pattern=separator_between_seq)
    n_separator_btw_seq_j_heavy <- stringi::stri_count_fixed(str=data[[j_call]], pattern=separator_between_seq)
    if (!all( n_separator_btw_seq_v_heavy == n_separator_btw_seq_j_heavy )) {
        stop("Requirement not met: one-to-one annotation-to-chain correspondence for both V and J (heavy)")
    }
    
    # check existence of additional columns, if specified
    if (!is.null(junc_len_heavy)) {
        if (!junc_len_heavy %in% colnames(data)) {
            stop("not all junc_len_heavy found as column(s) in data")
        }
    }
    if (!is.null(junc_len_light)) {
        if (!junc_len_light %in% colnames(data)) {
            stop("not all junc_len_light found as column(s) in data")
        }
    }
    
    # # STEP 0: Check V_CALL and J_CALL columns
    # chek na(s)
    int_nrow <- nrow(data)
    data <- data[!is.na(data[[v_call]]), ]
    data <- data[!is.na(data[[j_call]]), ]
    if (single_cell_with_light) {
        data <- data[!is.na(data[[v_call_light]]), ]
        data <- data[!is.na(data[[j_call_light]]), ]
    }
    if (!is.null(junc_len_heavy)) {
        data <- data[!is.na(data[[junc_len_heavy]]), ]
    }
    if (!is.null(junc_len_light)) {
        data <- data[!is.na(data[[junc_len_light]]), ]
    }
    fin_nrow <- nrow(data)
    if (int_nrow - fin_nrow > 0) {
        if (single_cell_with_light) {
            msg <- paste0("NA(s) found in one or more of { ", 
                          v_call, ", ", j_call, ", ", v_call_light, ", ", j_call_light, 
                          ifelse(is.null(junc_len_heavy), "", ", "), junc_len_heavy,
                          ifelse(is.null(junc_len_light), "", ", "), junc_len_light,
                          " } columns. ", int_nrow - fin_nrow, " sequence(s) removed.\n")
        } else {
            msg <- paste0("NA(s) found in one or more of { ", 
                          v_call, ", ", j_call, ifelse(is.null(junc_len_heavy), "", ", "), junc_len_heavy,
                          " } columns. ", int_nrow - fin_nrow, " sequence(s) removed.\n")
        }
        warning(msg)
    }
    
    ### expand
    
    # A1,A2;A3             [2 chains; 2/1 annotations per chain]
    # B1,B2;B1,B2;B3;B1,B3 [4 chains; 2/2/1/2 annotations per chain]
    
    # want 2*(2+2+1+2) + 1*(2+2+1+2) = (2+1)*(2+2+1+2) = 21
    
    # assumes that there is only one junction length per chain
    # (unlike annotations, which a chain could have multiple)
    
    # NULL will disappear when doing c()
    # c(NULL,NULL) gives NULL still
    cols_for_grouping_heavy <- c(v_call, j_call, junc_len_heavy)
    cols_for_grouping_light <- c(v_call_light, j_call_light, junc_len_light)
    
    # cols cannot be factor
    if (any( sapply(cols_for_grouping_heavy, function(x)class(data[[x]])) == "factor" )) {
        stop("one or more of { ", v_call, ", ", j_call,  
             ifelse(is.null(junc_len_heavy), " ", ", "), junc_len_heavy, 
             "} is factor. Must be character.\n")
    }
    if (!is.null(cols_for_grouping_light)) {
        if (any( sapply(cols_for_grouping_light, function(x)class(data[[x]])) == "factor" )) {
            stop("one or more of { ", v_call_light, ", ", j_call_light,  
                 ifelse(is.null(junc_len_light), " ", ", "), junc_len_light, 
                 "} is factor. Must be character.\n")
        }  
    }
    
    exp_lst <- vector(mode="list", length=nrow(data))
    
    for (i_orig in 1:nrow(data)) {
        
        #if (i_orig %% 1000 == 0) { cat(i_orig, "\n") }
        
        ## heavy
        cur_heavy_chains <- data[i_orig, cols_for_grouping_heavy]
        # matrix
        # each column is a chain
        # rows: V, J, (L)
        cur_heavy_chains_split <- stringi::stri_split_fixed(str=cur_heavy_chains, pattern=separator_between_seq, simplify=T)
        
        
        if (first) {
            # get gene name
            cur_heavy_chains_split[1:2, ] <- getGene(cur_heavy_chains_split[1:2, ], first=TRUE)
            cur_heavy_chains_split_total <- ncol(cur_heavy_chains_split)
        } else {
            # get gene name
            cur_heavy_chains_split[1:2, ] <- getGene(cur_heavy_chains_split[1:2, ], first=FALSE)
            # number of expansion per chain
            cur_heavy_chains_split_num <- sapply(1:ncol(cur_heavy_chains_split), function(i_chain){ 
                (stringi::stri_count_fixed(str=cur_heavy_chains_split[1,i_chain], pattern=separator_within_seq)+1) * 
                    (stringi::stri_count_fixed(str=cur_heavy_chains_split[2,i_chain], pattern=separator_within_seq)+1) })
            # number of expanded
            cur_heavy_chains_split_total <- sum(cur_heavy_chains_split_num)
        }
        
        cur_heavy_chains_exp <- rep(NA, length(cur_heavy_chains_split_total))
        i_exp <- 1
        for (i_chain in 1:ncol(cur_heavy_chains_split)) {
            v_i_chain <- stringi::stri_split_fixed(str=cur_heavy_chains_split[1, i_chain], pattern=separator_within_seq, simplify=F)[[1]]
            j_i_chain <- stringi::stri_split_fixed(str=cur_heavy_chains_split[2, i_chain], pattern=separator_within_seq, simplify=F)[[1]]
            
            for (ii_v in 1:length(v_i_chain)) {
                for (ii_j in 1:length(j_i_chain)) {
                    cur_heavy_chains_exp[i_exp] <- paste( v_i_chain[ii_v], j_i_chain[ii_j], 
                                                          ifelse(nrow(cur_heavy_chains_split)==3, cur_heavy_chains_split[3, i_chain], ""), 
                                                          sep="@")
                    i_exp <- i_exp+1
                }
            }
        }
        
        
        if (single_cell_with_light) {
            
            ## light
            
            cur_light_chains <- data[i_orig, cols_for_grouping_light]
            # matrix
            # each column is a chain
            # rows: V, J, (L)
            cur_light_chains_split <- stringi::stri_split_fixed(str=cur_light_chains, pattern=separator_between_seq, simplify=T)
            
            if (first) {
                # get gene name
                cur_light_chains_split[1:2, ] <- getGene(cur_light_chains_split[1:2, ], first=TRUE)
                cur_light_chains_split_total <- ncol(cur_light_chains_split)
            } else {
                # get gene name
                cur_light_chains_split[1:2, ] <- getGene(cur_light_chains_split[1:2, ], first=FALSE)
                # number of expansion per chain
                cur_light_chains_split_num <- sapply(1:ncol(cur_light_chains_split), function(i_chain){ 
                    (stringi::stri_count_fixed(str=cur_light_chains_split[1,i_chain], pattern=separator_within_seq)+1) * 
                        (stringi::stri_count_fixed(str=cur_light_chains_split[2,i_chain], pattern=separator_within_seq)+1) })
                # number of expanded
                cur_light_chains_split_total <- sum(cur_light_chains_split_num)
            }
            
            cur_light_chains_exp <- rep(NA, length(cur_light_chains_split_total))
            i_exp <- 1
            for (i_chain in 1:ncol(cur_light_chains_split)) {
                v_i_chain <- stringi::stri_split_fixed(str=cur_light_chains_split[1, i_chain], pattern=separator_within_seq, simplify=F)[[1]]
                j_i_chain <- stringi::stri_split_fixed(str=cur_light_chains_split[2, i_chain], pattern=separator_within_seq, simplify=F)[[1]]
                
                for (ii_v in 1:length(v_i_chain)) {
                    for (ii_j in 1:length(j_i_chain)) {
                        cur_light_chains_exp[i_exp] <- paste( v_i_chain[ii_v], j_i_chain[ii_j], 
                                                              ifelse(nrow(cur_light_chains_split)==3, cur_light_chains_split[3, i_chain], ""), 
                                                              sep="@")
                        i_exp <- i_exp+1
                    }
                }
            }
            
            ## pair with heavy
            cur_all <- expand.grid(cur_heavy_chains_exp, cur_light_chains_exp, stringsAsFactors=FALSE)
            cur_all_cat <- sapply(1:nrow(cur_all), function(x){paste0(cur_all[x,], collapse="@")})
        } else {
            cur_all_cat <- cur_heavy_chains_exp
        }
        
        exp_lst[[i_orig]] <- cur_all_cat
        
    }
    
    exp_uniq <- sort(unique(unlist(exp_lst)))
    n_cells_or_seqs <- length(exp_lst)
    
    # notes on implementation
    
    # regular/dense matrix is more straightforward to implement but very costly memory-wise
    # sparse matrix is less straightforward to implement but way more memory efficient
    
    # sparse matrix is very slow to modify to on-the-fly (using a loop like for dense matrix)
    # way faster to construct in one go
    
    # (DO NOT DELETE)
    # for illustrating the concept 
    # this is the way to go if using regular matrix (memory-intensive)
    # same concept implemented using sparse matrix
    
    # mtx_cell_VJL <- matrix(0, nrow=nrow(data), ncol=length(exp_uniq))
    # colnames(mtx_cell_VJL) <- exp_uniq
    # 
    # mtx_adj <- matrix(0, nrow=length(exp_uniq), ncol=length(exp_uniq))
    # rownames(mtx_adj) <- exp_uniq
    # colnames(mtx_adj) <- exp_uniq
    # 
    # for (i_cell in 1:length(exp_lst)) {
    #     #if (i_cell %% 1000 == 0) { cat(i_cell, "\n") }
    #     cur_uniq <- unique(exp_lst[[i_cell]])
    #     mtx_cell_VJL[i_cell, cur_uniq] <- 1
    #     mtx_adj[cur_uniq, cur_uniq] <- 1
    # }
    
    # actual implementation using sparse matrix from Matrix package
    
    ### matrix indicating relationship between cell and VJ(L) combinations
    # row: cell
    # col: unique heavy VJ(L) (and light VJ(L))
    
    # row indices
    m1_i <- lapply(1:length(exp_lst), function(i){rep(i, length(exp_lst[[i]]))})
    m1_i_v <- unlist(m1_i)
    
    # column indices
    m1_j <- lapply(exp_lst, function(x){
        idx <- match(unique(x), exp_uniq)
        #stopifnot( all.equal( exp_uniq[idx], unique(x) ) )
    })
    m1_j_v <- unlist(m1_j)
    
    stopifnot( length(m1_i_v) == length(m1_j_v) )
    
    # no particular need for this to be not of class "nsparseMatrix"
    # so no need to specify x=rep(1, length(m1_i))
    # not specifying makes it even more space-efficient
    mtx_cell_VJL <- Matrix::sparseMatrix(i=m1_i_v, j=m1_j_v, 
                                        dims=c(n_cells_or_seqs, length(exp_uniq)), 
                                        symmetric=F, triangular=F, index1=T, 
                                        dimnames=list(NULL, exp_uniq))
    
    ### adjacency matrix
    # row and col: unique heavy VJ(L) (and light VJ(L))
    
    # row indices
    m2_i <- lapply(m1_j, function(x){ rep(x, each=length(x)) })
    m2_i_v <- unlist(m2_i)
    
    # col indices
    m2_j <- lapply(m1_j, function(x){ rep(x, times=length(x)) })
    m2_j_v <- unlist(m2_j)
    
    stopifnot( length(m2_i_v) == length(m2_j_v) )
    
    # important: x must be specified for mtx_adj in order to make it not of class "nsparseMatrix"
    # this is because igraph accepts sparse matrix from Matrix but not the "pattern" matrices variant
    mtx_adj <- Matrix::sparseMatrix(i=m2_i_v, j=m2_j_v, x=rep(1,length(m2_i_v)), 
                                   dims=c(length(exp_uniq), length(exp_uniq)), 
                                   symmetric=F, triangular=F, index1=T, 
                                   dimnames=list(exp_uniq, exp_uniq))
    mtx_adj <- mtx_adj>0
    
    rm(m1_i, m1_j, m2_i, m2_j, m1_i_v, m1_j_v, m2_i_v, m2_j_v, exp_lst)
    
    ### identify connected components based on adjcencey matrix
    # this is the grouping
    # source: https://stackoverflow.com/questions/35772846/obtaining-connected-components-in-r
    
    g <- igraph::graph_from_adjacency_matrix(adjmatrix=mtx_adj, mode="undirected", diag=FALSE)
    #plot(g, vertex.size=10, vertex.label.cex=1, vertex.color="skyblue", vertex.label.color="black", vertex.frame.color="transparent", edge.arrow.mode=0)
    
    connected <- igraph::components(g)
    VJL_groups <- igraph::groups(connected)
    names(VJL_groups) <- paste0("G", 1:length(VJL_groups))
    
    ### identify cells associated with each connected component (grouping)
    
    # each entry corresponds to a group/partition
    # each element within an entry is a cell
    
    cellIdx_byGroup_lst <- lapply(VJL_groups, function(x){ 
        if (length(x)>1) {
            # matrix
            # important to specify rowSums from Matrix package
            # base::rowSums will NOT work
            cell_idx <- which(Matrix::rowSums(mtx_cell_VJL[, x])>0)
        } else {
            # vector
            cell_idx <- which(mtx_cell_VJL[, x]>0)
        }
        return(cell_idx)
    })
    
    # sanity check: there should be perfect/disjoint partitioning 
    # (each cell has exactly one group assignment)
    stopifnot( n_cells_or_seqs == length(unique(unlist(cellIdx_byGroup_lst))) )
    
    # assign
    data$VJ_GROUP <- NA
    for (i in 1:length(cellIdx_byGroup_lst)) {
        data[["VJ_GROUP"]][cellIdx_byGroup_lst[[i]]] <- names(VJL_groups)[i]
    }
    stopifnot(!any(is.na(data[["VJ_GROUP"]])))
    
    return(data)
}


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
