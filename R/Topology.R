# Ig lineage topology analysis

#' @include Classes.R
NULL


#### Graph analysis functions ####

#' Generate subtree summary statistics for a tree
#'
#' \code{summarizeSubtrees} calculates summary statistics for each node of a tree. Includes
#' both node properties and subtree properties.
#'
#' @param    graph   igraph object containing an annotated lineage tree.
#' @param    fields  annotation fields to add to the output.
#' @param    root    name of the root (germline) node.
#' 
#' @return   A data.frame with columns: 
#'           \itemize{
#'             \item  \code{NAME}:             node name.
#'             \item  \code{PARENT}:           name of the parent node.
#'             \item  \code{OUTDEGREE}:        number of edges leading from the node.
#'             \item  \code{SIZE}:             total number of nodes within the subtree rooted 
#'                                             at the node.
#'             \item  \code{DEPTH}:            the depth of the subtree that is rooted at 
#'                                             the node.
#'             \item  \code{PATHLENGTH}:       the maximum pathlength beneath the node.
#'             \item  \code{OUTDEGREE_NORM}:   \code{OUTDEGREE} normalized by the total 
#'                                             number of edges.
#'             \item  \code{SIZE_NORM}:        \code{SIZE} normalized by the largest
#'                                             subtree size (the germline).
#'             \item  \code{DEPTH_NORM}:       \code{DEPTH} normalized by the largest
#'                                             subtree depth (the germline).
#'             \item  \code{PATHLENGTH_NORM}:  \code{PATHLEGNTH} normalized by the largest
#'                                             subtree pathlength (the germline).
#'           }
#'           An additional column corresponding to the value of \code{field} is added when
#'           specified.
#' 
#' @seealso  See \link{buildPhylipLineage} for generating input trees. 
#'           See \link{getPathLengths} for calculating path length to nodes.
#' 
#' @examples
#' # Define and plot example graph
#' library(igraph)
#' graph <- ExampleTrees[[23]]
#' plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE)
#' 
#' # Summarize tree
#' summarizeSubtrees(graph, fields="ISOTYPE", root="Germline")
#' 
#' @export
summarizeSubtrees <- function(graph, fields=NULL, root="Germline") {
    ## DEBUG
    # root="Germline"; fields=NULL
    # TODO:  should probably include a means to exclude inferred from substree size
    
    # Define node attribute data.frame    
    node_df <- data.frame(NAME=V(graph)$name, stringsAsFactors=F)
    for (f in fields) { 
        node_df[[f]] <- vertex_attr(graph, name=f) 
    }

    # Get edges
    edges <- igraph::as_edgelist(graph)
    # Get unweighted paths
    paths_step <- suppressWarnings(igraph::distances(graph, mode="out", algorithm="unweighted"))
    paths_step[!is.finite(paths_step)] <- NA
    # Get weighted paths
    paths_length <- igraph::distances(graph, mode="out", algorithm="dijkstra")
    paths_length[!is.finite(paths_length)] <- NA
    
    # Define each node's parent
    node_df$PARENT <- edges[, 1][match(node_df$NAME, edges[, 2])]
    # Define each node's outdegree
    node_df$OUTDEGREE <- igraph::degree(graph, mode="out")
    # Define the number of nodes in each subtree (child count + 1)
    node_df$SIZE <- apply(paths_step, 1, function(x) length(na.omit(x)))
    # Define number of levels below each node
    node_df$DEPTH <- apply(paths_step, 1, max, na.rm=TRUE) + 1
    # Define the maximum shortest path length (genetic distance) to a leaf from each node
    node_df$PATHLENGTH <- apply(paths_length, 1, max, na.rm=TRUE)
    
    # Normalize
    node_df <- node_df %>%
        dplyr::mutate_(OUTDEGREE_NORM=interp(~x/sum(x, na.rm=TRUE), x=as.name("OUTDEGREE")),
                       SIZE_NORM=interp(~x/max(x, na.rm=TRUE), x=as.name("SIZE")),
                       DEPTH_NORM=interp(~x/max(x, na.rm=TRUE), x=as.name("DEPTH")),
                       PATHLENGTH_NORM=interp(~x/max(x, na.rm=TRUE), x=as.name("PATHLENGTH")))

    return(node_df)
}


#' Calculate path lengths from the tree root
#'
#' \code{getPathLengths} calculates the unweighted (number of steps) and weighted (distance) 
#' path lengths from the root of a lineage tree.
#'
#' @param    graph     igraph object containing an annotated lineage tree.
#' @param    root      name of the root (germline) node.
#' @param    field     annotation field to use for exclusion of nodes from step count.
#' @param    exclude   annotation values specifying which nodes to exclude from step count. 
#'                     If \code{NULL} consider all nodes. This does not affect the weighted
#'                     (distance) path length calculation.
#'                     
#' @return   A data.frame with columns:
#'           \itemize{
#'             \item  \code{NAME}:      node name
#'             \item  \code{STEPS}:     path length as the number of nodes traversed
#'             \item  \code{DISTANCE}:  path length as the sum of edge weights
#'           }
#' 
#' @seealso  See \link{buildPhylipLineage} for generating input trees. 
#' 
#' @examples
#' # Define and plot example graph
#' library(igraph)
#' graph <- ExampleTrees[[23]]
#' plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE)
#' 
#' # Consider all nodes
#' getPathLengths(graph, root="Germline")
#' 
#' # Exclude nodes without an isotype annotation from step count
#' getPathLengths(graph, root="Germline", field="ISOTYPE", exclude=NA)
#' 
#' @export
getPathLengths <- function(graph, root="Germline", field=NULL, exclude=NULL) {
    # Define path length data.frame
    path_df <- data.frame(NAME=V(graph)$name, stringsAsFactors=FALSE)

    # Get indices of excluded vertices
    skip_idx <- which(path_df$NAME == root)
    if (!is.null(field)) {
        g <- vertex_attr(graph, name=field)
        skip_idx <- union(skip_idx, which(g %in% exclude))
    }
    
    # Get paths
    step_list <- shortest_paths(graph, root, mode="out", weights=NA, output="vpath")
    step_list <- step_list$vpath
    
    # Get path lengths
    for (i in 1:length(step_list)) {
        v <- step_list[[i]]
        path_df[i, "STEPS"] <- sum(!(v %in% skip_idx)) 
        path_df[i, "DISTANCE"] <- sum(E(graph, path=v)$weight)
    }
    
    return(path_df)
}


#' Retrieve the first non-root node of a lineage tree
#' 
#' \code{getMRCA} returns the set of lineage tree nodes with the minimum weighted or 
#' unweighted path length from the root (germline) of the lineage tree, allowing for 
#' exclusion of specific groups of nodes.
#'
#' @param    graph    igraph object containing an annotated lineage tree.
#' @param    path     string defining whether to use unweighted (steps) or weighted (distance) 
#'                    measures for determining the founder node set.. 
#' @param    root     name of the root (germline) node.
#' @param    field    annotation field to use for both unweighted path length exclusion and
#'                    consideration as an MRCA node. If \code{NULL} do not exclude any nodes.
#' @param    exclude  vector of annotation values in \code{field} to exclude from the potential 
#'                    MRCA set. If \code{NULL} do not exclude any nodes. Has no effect if 
#'                    \code{field=NULL}.
#'                    
#' @return   A data.frame of the MRCA node(s) containing the columns:
#'           \itemize{
#'             \item  \code{NAME}:      node name
#'             \item  \code{STEPS}:     path length as the number of nodes traversed
#'             \item  \code{DISTANCE}:  path length as the sum of edge weights
#'           }
#'           Along with additional columns corresponding to the 
#'           annotations of the input graph.
#'           
#' @seealso  Path lengths are determined with \link{getPathLengths}.
#' 
#' @examples
#' # Define and plot example graph
#' library(igraph)
#' graph <- ExampleTrees[[23]]
#' plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE)
#' 
#' # Use unweighted path length and do not exclude any nodes
#' getMRCA(graph, path="steps", root="Germline")
#'
#' # Exclude nodes without an isotype annotation and use weighted path length
#' getMRCA(graph, path="distance", root="Germline", field="ISOTYPE", exclude=NA)
#' 
#' @export
getMRCA <- function(graph, path=c("distance", "steps"), root="Germline", 
                    field=NULL, exclude=NULL) {
    # Check arguments
    path <- match.arg(path)
    
    # Get distance from root
    path_df <- getPathLengths(graph, root=root, field=field, exclude=exclude)
    
    # Get indices of excluded vertices
    skip_idx <- which(path_df$NAME == root)
    if (!is.null(field)) {
        g <- vertex_attr(graph, name=field)
        skip_idx <- union(skip_idx, which(g %in% exclude))
    }
    
    # Get founder nodes
    if (path == "distance") { 
        path_len <- setNames(path_df$DISTANCE, 1:nrow(path_df))
    } else if (path == "steps") {
        path_len <- setNames(path_df$STEPS, 1:nrow(path_df))
    } else {
        stop("Invalid value for 'path' parameter. Must be one of c('distance', 'steps').\n")
    }
    
    path_len <- path_len[-skip_idx]
    root_idx <- as.numeric(names(path_len)[which(path_len == min(path_len))])
    root_df <- igraph::as_data_frame(graph, what="vertices")[root_idx, ]
    root_df$STEPS <- path_df$STEPS[root_idx]
    root_df$DISTANCE <- path_df$DISTANCE[root_idx]
    
    # Switch name column to uppercase
    names(root_df)[names(root_df) == "name"] <- "NAME"
    
    return(root_df)
}


#' Tabulate the number of edges between annotations within a lineage tree
#'
#' \code{tableEdges} creates a table of the total number of connections (edges) for each 
#' unique pair of annotations within a tree over all nodes.
#' 
#' @param    graph     igraph object containing an annotated lineage tree.
#' @param    field     string defining the annotation field to count.
#' @param    indirect  if \code{FALSE} count direct connections (edges) only. If 
#'                     \code{TRUE} walk through any nodes with annotations specified in 
#'                     the \code{argument} to count indirect connections. Specifying
#'                     \code{indirect=TRUE} with \code{exclude=NULL} will have no effect.
#' @param    exclude   vector of strings defining \code{field} values to exclude from counts.
#'                     Edges that either start or end with the specified annotations will not
#'                     be counted. If \code{NULL} count all edges.
#'                     
#' @return   A data.frame defining total annotation connections in the tree with columns:
#'           \itemize{
#'             \item  \code{PARENT}:  parent annotation
#'             \item  \code{CHILD}:   child annotation
#'             \item  \code{COUNT}:   count of edges for the parent-child relationship
#'           }
#'           
#' @seealso  See \link{testEdges} for performed a permutation test on edge relationships.
#'           
#' @examples
#' # Define and plot example graph
#' library(igraph)
#' graph <- ExampleTrees[[23]]
#' plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE)
#' 
#' # Count direct edges between isotypes including inferred nodes
#' tableEdges(graph, "ISOTYPE", exclude=NULL)
#' 
#' # Count direct edges excluding edges to and from inferred nodes
#' tableEdges(graph, "ISOTYPE", exclude=NA)
#' 
#' # Count indirect edges walking through inferred nodes
#' tableEdges(graph, "ISOTYPE", indirect=TRUE, exclude=NA)
#' 
#' @export
tableEdges <- function(graph, field, indirect=FALSE, exclude=c("Germline", NA)) {
    # Function to retrieve the name if x is exactly one vertex index and NULL otherwise
    .getSingleVertex <- function(x) {
        if (length(x) == 1) { 
            vertex_attr(graph, name=field, index=x[1]) 
        } else { 
            NULL 
        }
    }

    if (indirect) {
        # Get indices of excluded and retained vertices
        if (!is.null(exclude)) {
            f <- vertex_attr(graph, name=field)
            skip_idx <- which(f %in% exclude)
            keep_idx <- as.numeric(V(graph))[-skip_idx]
        } else {
            skip_idx <- NULL
            keep_idx <- as.numeric(V(graph))
        }
        
        # Iterate over nodes and count indirect parent-child connections
        edge_list <- list()
        for (i in keep_idx) { 
            # Get parent annotation
            parent <- vertex_attr(graph, name=field, index=i)
            
            # Get indirect child node annotations
            step_list <- suppressWarnings(shortest_paths(graph, V(graph)[i], mode="out", weights=NA, output="vpath"))
            step_list <- unique(lapply(step_list$vpath, function(x) x[!(x %in% c(i, skip_idx))]))
            children <- unlist(lapply(step_list, .getSingleVertex))
                
            # Define data.frame of connections
            if (length(children) > 0) {
                edge_list[[i]] <- data.frame("PARENT"=parent, "CHILD"=children, 
                                             stringsAsFactors=FALSE)
            }
        }
        
        # Merge edge list into data.frame
        edge_df <- bind_rows(edge_list)        
    }
    else {
        # Get adjacency list
        edge_mat <- as_edgelist(graph, names=FALSE)
        edge_mat <- vertex_attr(graph, name=field, index=edge_mat)
        edge_mat <- matrix(edge_mat, ncol=2, dimnames=list(NULL, c("PARENT", "CHILD")))

        # Build and subset edge data.frame
        edge_df <- as.data.frame(edge_mat, stringsAsFactors=FALSE)
        edge_df <- edge_df[!(edge_df$PARENT %in% exclude) & !(edge_df$CHILD %in% exclude), ]
    }
    
    # Count edges
    edge_tab <- edge_df %>%
        group_by_("PARENT", "CHILD") %>%
        dplyr::summarize(COUNT=n())

    return(edge_tab)
}


#' Permute the node labels of a tree
#' 
#' \code{permuteLabels} permutes the node annotations of a lineage tree.
#'
#' @param    graph    igraph object containing an annotated lineage tree.
#' @param    field    string defining the annotation field to permute.
#' @param    exclude  vector of strings defining \code{field} values to exclude 
#'                    from permutation.
#' 
#' @return   A modified igraph object with vertex annotations permuted.
#' 
#' @seealso  \link{testEdges}.
#' 
#' @examples
#' # Define and plot example graph
#' library(igraph)
#' graph <- ExampleTrees[[23]]
#' plot(graph, layout=layout_as_tree, vertex.label=V(graph)$ISOTYPE)
#' 
#' # Permute annotations and plot new tree
#' g <- permuteLabels(graph, "ISOTYPE")
#' plot(g, layout=layout_as_tree, vertex.label=V(g)$ISOTYPE)
#' 
#' @export
permuteLabels <- function(graph, field, exclude=c("Germline", NA)) {
    # Determine which nodes to permute
    labels <- vertex_attr(graph, name=field)
    i <- which(!(labels %in% exclude))
    
    # Return input on insufficient number of nodes
    if (length(i) < 2) { 
        warning("Only 1 node to permute\n")
        return(graph) 
    }
    
    # Sample and reassign field values
    s <- sample(i)
    perm <- set_vertex_attr(graph, name=field, index=i, value=labels[s])
    
    return(perm)
}



#### Test functions ####


#' Tests for MRCA annotation enrichment in lineage trees
#' 
#' \code{testMRCA} performs a permutation test on a set of lineage trees to determine
#' the significance of an annotation's association with the MRCA position of the lineage
#' trees.
#' 
#' @param    graphs    list of igraph object containing annotated lineage trees.
#' @param    field     string defining the annotation field to test.
#' @param    root      name of the root (germline) node.
#' @param    exclude   vector of strings defining \code{field} values to exclude from the
#'                     set of potential founder annotations.
#' @param    nperm     number of permutations to perform.
#' @param    progress  if \code{TRUE} show a progress bar.
#' 
#' @return   An \link{MRCATest} object containing the test results and permutation
#'           realizations.
#'           
#' @seealso  Uses \link{getMRCA} and \link{getPathLengths}. 
#'           See \link{plotMRCATest} for plotting the permutation distributions.
#'           
#' @examples
#' # Perform MRCA test on isotypes
#' x <- testMRCA(ExampleTrees, "ISOTYPE", nperm=100)
#' print(x)
#' 
#' @export
testMRCA <- function(graphs, field, root="Germline", exclude=c("Germline", NA), 
                     nperm=200, progress=FALSE) {
    # Function to resolve ambiguous founders
    # @param  x      data.frame from getMRCA
    # @param  field  annotation field
    .resolveMRCA <- function(x, field) {
        x %>% filter_(interp(~!duplicated(y), y=as.name(field))) %>%
            filter_(interp(~length(y) == 1, y=as.name(field)))
    }
    
    # Function to count MRCAs
    # @param  x        list of graphs
    # @param  field    annotation field
    # @param  exclude  vector of annotation values to exclude
    .countMRCA <- function(x, field, exclude) {
        # Get MRCAs
        mrca_list <- lapply(x, getMRCA, path="distance", field=field, 
                            exclude=exclude)
        # Resolve ambiguous MRCAs
        mrca_list <- lapply(mrca_list, .resolveMRCA, field=field)
        # Summarize MRCA counts
        mrca_sum <- bind_rows(mrca_list, .id="GRAPH") %>%
            select_("GRAPH", field) %>%
            rename_("ANNOTATION"=field) %>%
            group_by_("ANNOTATION") %>%
            dplyr::summarize(COUNT=n())
        
        return(mrca_sum)
    }
    
    # Assign numeric names if graphs is an unnamed list
    if (is.null(names(graphs))) { names(graphs) <- 1:length(graphs) }
    
    # Summarize observed MRCA counts
    obs_sum <- .countMRCA(graphs, field=field, exclude=exclude)

    # Generate edge null distribution via permutation
    if (progress) { 
        cat("-> PERMUTING TREES\n")
        pb <- txtProgressBar(min=0, max=nperm, initial=0, width=40, style=3) 
    }
    perm_list <- list()
    for (i in 1:nperm) {
        # Permute labels
        tmp_list <- lapply(graphs, permuteLabels, field=field, exclude=exclude)
        # Summarize MRCA counts
        tmp_sum <- .countMRCA(tmp_list, field=field, exclude=exclude)
        # Update permutation set
        tmp_sum$ITER <- i
        perm_list[[i]] <- tmp_sum
        
        if (progress) { setTxtProgressBar(pb, i) }
    }
    cat("\n")
    perm_sum <- bind_rows(perm_list)
    
    # Test observed against permutation distribution
    for (i in 1:nrow(obs_sum)) {
        x <- obs_sum$ANNOTATION[i]
        # Annotation count distribution
        d <- perm_sum$COUNT[perm_sum$ANNOTATION == x]
        # Expected mean
        obs_sum[i, "EXPECTED"] <- mean(d)
        # P-value for observed > expected
        f <- ecdf(d)
        obs_sum[i, "PVALUE"] <- 1 - f(obs_sum$COUNT[i])
    }
 
    # Generate return object
    mrca_test <- new("MRCATest", 
                     tests=as.data.frame(obs_sum), 
                     permutations=as.data.frame(perm_sum),
                     nperm=nperm)
    
    return(mrca_test)
}


#' Tests for parent-child annotation enchrichment in lineage trees
#' 
#' \code{testEdges} performs a permutation test on a set of lineage trees to determine
#' the significance of an annotation's association with parent-child relationships.
#'
#' @param    graphs    list of igraph objects with vertex annotations.
#' @param    field     string defining the annotation field to permute.
#' @param    exclude   vector of strings defining \code{field} values to exclude from 
#'                     permutation.
#' @param    nperm     number of permutations to perform.
#' @param    progress  if \code{TRUE} show a progress bar.
#' 
#' @return   An \link{EdgeTest} object containing the test results and permutation
#'           realizations.
#' 
#' @seealso  Uses \link{tableEdges} and \link{permuteLabels}. 
#'           See \link{plotEdgeTest} for plotting the permutation distributions.
#'           
#' @examples
#' # Perform edge test on isotypes
#' x <- testEdges(ExampleTrees, "ISOTYPE", nperm=100)
#' print(x)
#' 
#' @export
testEdges <- function(graphs, field, exclude=c("Germline", NA), nperm=200, 
                      progress=FALSE) {
    ## DEBUG
    # field="isotype"; exclude=c("Germline", NA); nperm=200
    
    # Assign numeric names if graphs is an unnamed list
    if (is.null(names(graphs))) { names(graphs) <- 1:length(graphs) }
    
    # Function to count edge annotations
    # @param  x        list of graphs
    # @param  field    annotation field
    # @param  exclude  vector of annotation values to exclude
    .countEdges <- function(x, field, exclude) {
        edge_list <- lapply(x, tableEdges, field=field, exclude=exclude)
        edge_sum <- bind_rows(edge_list) %>%
            group_by_("PARENT", "CHILD") %>%
            dplyr::summarize_(COUNT=interp(~sum(x, na.rm=TRUE), x=as.name("COUNT")))
        return(edge_sum)
    }
    
    # Count edges of observed data
    obs_sum <- .countEdges(graphs, field, exclude)

    # Generate edge null distribution via permutation
    if (progress) { 
        cat("-> PERMUTING TREES\n")
        pb <- txtProgressBar(min=0, max=nperm, initial=0, width=40, style=3) 
    }
    perm_list <- list()
    for (i in 1:nperm) {
        # Permute annotations
        tmp_list <- lapply(graphs, permuteLabels, field=field, exclude=exclude)
        # Count edges
        tmp_sum <- .countEdges(tmp_list, field, exclude)
        # Update permutation set
        tmp_sum$ITER <- i
        perm_list[[i]] <- tmp_sum
        
        if (progress) { setTxtProgressBar(pb, i) }
    }
    cat("\n")
    perm_sum <- bind_rows(perm_list)
    
    # Test observed against permutation distribution
    for (i in 1:nrow(obs_sum)) {
        x <- obs_sum$PARENT[i]
        y <- obs_sum$CHILD[i]
        # Edge count distribution
        d <- perm_sum$COUNT[perm_sum$PARENT == x & perm_sum$CHILD == y]
        # Expected mean
        obs_sum[i, "EXPECTED"] <- mean(d)
        # P-value for observed > expected
        f <- ecdf(d)
        obs_sum[i, "PVALUE"] <- 1 - f(obs_sum$COUNT[i])
    }
    
    # Generate return object
    edge_test <- new("EdgeTest", 
                 tests=as.data.frame(obs_sum), 
                 permutations=as.data.frame(perm_sum),
                 nperm=nperm)
    
    return(edge_test)
}


#### Plotting functions #####

#' Plot the results of an edge permutation test
#' 
#' \code{plotEdgeTest} plots the results of an edge permutation test performed with 
#' \code{testEdges} as either a histogram or cumulative distribution function.
#'
#' @param    data        \link{EdgeTest} object returned by \link{testEdges}.
#' @param    color       color of the histogram or lines.
#' @param    main_title  string specifying the plot title.
#' @param    style       type of plot to draw. One of:
#'                       \itemize{
#'                         \item \code{"histogram"}:  histogram of the edge count 
#'                                                    distribution with a red dotted line
#'                                                    denoting the observed value.
#'                         \item \code{"cdf"}:        cumulative distribution function 
#'                                                    of edge counts with a red dotted 
#'                                                    line denoting the observed value and
#'                                                    a blue dotted line indicating the 
#'                                                    p-value.
#'                       }
#' @param    silent      if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                       object; if \code{FALSE} draw the plot.
#' @param    ...         additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  See \link{testEdges} for performing the test.
#' 
#' @examples
#' # Perform edge test on isotypes
#' x <- testEdges(ExampleTrees, "ISOTYPE", nperm=100)
#' 
#' # Plot
#' plotEdgeTest(x, color="steelblue", style="hist")
#' plotEdgeTest(x, style="cdf")
#' 
#' @export
plotEdgeTest <- function(data, color="black", main_title="Edge Test", 
                         style=c("histogram", "cdf"), silent=FALSE, ...) {
    # Check arguments
    style <- match.arg(style)
    
    # Extract plot data
    obs_sum <- rename_(data@tests, "Parent"="PARENT", "Child"="CHILD")
    perm_sum <- rename_(data@permutations, "Parent"="PARENT", "Child"="CHILD")

    if (style == "histogram") {
        # Plot edge null distribution
        p1 <- ggplot(perm_sum, aes_string(x="COUNT")) +
            getBaseTheme() + 
            ggtitle(main_title) +
            xlab("Number of edges") +
            ylab("Number of realizations") + 
            geom_histogram(binwidth=1,
                           fill=color, color=NA) +
            geom_vline(data=obs_sum, aes_string(xintercept="COUNT"), 
                       color="firebrick", linetype=3, size=1.25) + 
            facet_grid("Child ~ Parent", labeller=label_both, scales="free")
    } else if (style == "cdf") {    
        # Plot ECDF of edge null distribution
        p1 <- ggplot(perm_sum, aes_string(x="COUNT")) +
            getBaseTheme() + 
            ggtitle(main_title) +
            xlab("Number of edges") +
            ylab("P-value") +
            stat_ecdf(color=color, size=1) +
            geom_vline(data=obs_sum, aes_string(xintercept="COUNT"), color="firebrick", 
                       linetype=3, size=1.25) + 
            geom_hline(data=obs_sum, aes_string(yintercept="PVALUE"), color="steelblue", 
                       linetype=3, size=1.25) + 
            facet_grid("Child ~ Parent", labeller=label_both, scales="free")
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}


#' Plot the results of a founder permutation test
#' 
#' \code{plotMRCATest} plots the results of a founder permutation test performed with 
#' \code{testMRCA}.
#'
#' @param    data        \link{MRCATest} object returned by \link{testMRCA}.
#' @param    color       color of the histogram or lines.
#' @param    main_title  string specifying the plot title.
#' @param    style       type of plot to draw. One of:
#'                       \itemize{
#'                         \item \code{"histogram"}:  histogram of the annotation count 
#'                                                    distribution with a red dotted line
#'                                                    denoting the observed value.
#'                         \item \code{"cdf"}:        cumulative distribution function 
#'                                                    of annotation counts with a red dotted 
#'                                                    line denoting the observed value and
#'                                                    a blue dotted line indicating the 
#'                                                    p-value.
#'                       }
#' @param    silent      if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                       object; if \code{FALSE} draw the plot.
#' @param    ...         additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  See \link{testEdges} for performing the test.
#' 
#' @examples
#' # Perform MRCA test on isotypes
#' x <- testMRCA(ExampleTrees, "ISOTYPE", nperm=100)
#' 
#' # Plot
#' plotMRCATest(x, color="steelblue", style="hist")
#' plotMRCATest(x, style="cdf")
#' 
#' @export
plotMRCATest <- function(data, color="black", main_title="MRCA Test", 
                         style=c("histogram", "cdf"), silent=FALSE, ...) {
    # Check arguments
    style <- match.arg(style)
    
    # Extract plot data
    obs_sum <- rename_(data@tests, "Annotation"="ANNOTATION")
    perm_sum <- rename_(data@permutations, "Annotation"="ANNOTATION")
    
    if (style == "histogram") {
        # Plot MRCA null distribution
        p1 <- ggplot(perm_sum, aes_string(x="COUNT")) +
            getBaseTheme() + 
            ggtitle(main_title) +
            xlab("Number of MRCAs") +
            ylab("Number of realizations") + 
            geom_histogram(fill=color, color=NA, binwidth=1) +
            geom_vline(data=obs_sum, aes_string(xintercept="COUNT"), 
                       color="firebrick", linetype=3, size=1.25) + 
            facet_wrap("Annotation", ncol=1, scales="free_y")
    } else if (style == "cdf") {
        # Plot ECDF of MRCA null distribution
        p1 <- ggplot(perm_sum, aes_string(x="COUNT")) +
            getBaseTheme() + 
            ggtitle(main_title) +
            xlab("Number of MRCAs") +
            ylab("P-value") +
            stat_ecdf(color=color, size=1) +
            geom_vline(data=obs_sum, aes_string(xintercept="COUNT"), 
                       color="firebrick", linetype=3, size=1.25) + 
            geom_hline(data=obs_sum, aes_string(yintercept="PVALUE"), 
                       color="steelblue", linetype=3, size=1.25) + 
            facet_wrap("Annotation", nrow=1, scales="free_y")
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}


#' Plots subtree statistics for multiple trees
#' 
#' \code{plotSubtree} plots distributions of normalized subtree statistics for a 
#' set of lineage trees, broken down by annotation value.
#'
#' @param    graphs        list of igraph objects containing annotated lineage trees.
#' @param    field         string defining the annotation field.
#' @param    stat          string defining the subtree statistic to plot. One of:
#'                         \itemize{
#'                           \item  \code{outdegree}:   distribution of normalized node 
#'                                                      outdegrees.
#'                           \item  \code{size}:        distribution of normalized subtree sizes.
#'                           \item  \code{depth}:       distribution of subtree depths.
#'                           \item  \code{pathlength}:  distribution of maximum pathlength 
#'                                                      beneath nodes.
#'                         }
#' @param    root          name of the root (germline) node.
#' @param    exclude       vector of strings defining \code{field} values to exclude from
#'                         plotting.
#' @param    colors        named vector of colors for values in \code{field}, with 
#'                         names defining annotation names \code{field} column and values
#'                         being colors. Also controls the order in which values appear on the
#'                         plot. If \code{NULL} alphabetical ordering and a default color palette 
#'                         will be used.
#' @param    main_title    string specifying the plot title.
#' @param    legend_title  string specifying the legend title.
#' @param    style         string specifying the style of plot to draw. One of:
#'                         \itemize{
#'                           \item \code{"histogram"}:  histogram of the annotation count 
#'                                                      distribution with a red dotted line
#'                                                      denoting the observed value.
#'                           \item \code{"cdf"}:        cumulative distribution function 
#'                                                      of annotation counts with a red 
#'                                                      dotted line denoting the observed 
#'                                                      value and a blue dotted line 
#'                                                      indicating the p-value.
#'                         }
#' @param    silent        if \code{TRUE} do not draw the plot and just return the ggplot2 
#'                         object; if \code{FALSE} draw the plot.
#' @param    ...           additional arguments to pass to ggplot2::theme.
#'
#' @return   A \code{ggplot} object defining the plot.
#' 
#' @seealso  Subtree statistics are calculated with \link{summarizeSubtrees}.
#' 
#' @examples
#' # Plot boxplot of outdegree by sample
#' plotSubtrees(ExampleTrees, "SAMPLE", "out", main_title="Node outdegree", 
#'              style="b")
#'
#' # Plot boxplot of subtree by sample
#' plotSubtrees(ExampleTrees, "SAMPLE", "size", style="b")
#' 
#' # Plot violins of pathlength by isotype
#' plotSubtrees(ExampleTrees,  "ISOTYPE", "path", colors=IG_COLORS, 
#'              legend_title="Isotype", style="v")
#' 
#' # Plot violins of depth by isotype
#' plotSubtrees(ExampleTrees,  "ISOTYPE", "depth", style="v")
#' 
#' @export
plotSubtrees <- function(graphs, field, stat, root="Germline", exclude=c("Germline", NA), 
                         colors=NULL, main_title="Subtrees", legend_title="Annotation", 
                         style=c("box", "violin"), silent=FALSE, ...) {
    # Hack for visibility of special ggplot variables
    ..x.. <- ..y.. <- NULL
    
    ## DEBUG
    # graphs=ExampleTrees; field="ISOTYPE"; colors=IG_COLORS; main_title="Outdegree"; root="Germline"; exclude=c("Germline", NA); style="box"
    # Check arguments
    style <- match.arg(style, several.ok=FALSE)
    stat <- match.arg(stat, choices=c("outdegree", "size", "depth", "pathlength"), 
                      several.ok=FALSE)

    # Set stat column and axis labels
    if (stat == "outdegree") {
        stat_col <- "OUTDEGREE_NORM"
        y_lab <- "Node outdegree (percent of tree edges)"
    } else if (stat == "size") {
        stat_col <- "SIZE_NORM"
        y_lab <- "Substree size (percent of tree size)"
    } else if (stat == "depth") {
        stat_col <- "DEPTH_NORM"
        y_lab <- "Depth under node (percent of tree depth)"
    } else if (stat == "pathlength") {
        stat_col <- "PATHLENGTH_NORM"
        y_lab <- "Path length under node (percent of longest path)"
    } else {
        stop("Invalid value for 'stat'. How did you get here?")
    }
    
    # Assign numeric names if graphs is an unnamed list
    if (is.null(names(graphs))) { names(graphs) <- 1:length(graphs) }
    
    # Get subtree summarizes and filter excluded annotations
    sum_list <- lapply(graphs, summarizeSubtrees, fields=field, root=root)
    sum_df <- bind_rows(sum_list, .id="GRAPH") %>%
        filter_(interp(~!(x %in% exclude), x=as.name(field)),
                interp(~is.finite(x), x=as.name(stat_col)))
    
    # Set ordering based on color names
    if (!is.null(colors)) {
        # Assign missing levels to grey
        x <- unique(sum_df[[field]])
        x <- sort(x[!(x %in% names(colors))])
        if (length(x) > 0) {
            warning("The following annotations are missing from the 'colors' argument and will be assigned to grey: ", 
                    paste(x, collapse=" "))
            x <- setNames(rep("grey", length(x)), x)
            colors <- c(colors, x)
        }
        # Cast to factor
        sum_df[[field]] <- factor(sum_df[[field]], levels=names(colors))
    } else {
        sum_df[[field]] <- factor(sum_df[[field]])
    }
    
    # Make plot object
    p1 <- ggplot(sum_df, aes_string(x=field, y=stat_col)) + 
        getBaseTheme() + 
        ggtitle(main_title) + 
        xlab("") +
        ylab(y_lab) +
        scale_y_continuous(labels=percent) +
        expand_limits(y=0)
    
    # Add distributions style
    if (style == "box") {
        p1 <- p1 + geom_boxplot(aes_string(fill=field), width=0.7, alpha=0.8)
    } else if (style == "violin") {
        p1 <- p1 + geom_violin(aes_string(fill=field), adjust=1.5, scale="width", trim=T, 
                               width=0.7, alpha=0.8) +
            geom_errorbarh(aes(xmin=..x.. - 0.4, xmax=..x.. + 0.4), color="black", 
                           stat="summary", fun.y="mean", size=1.25, height=0, alpha=0.9)
    }

    # Set colors and legend
    if (!is.null(colors)) {
        p1 <- p1 + scale_fill_manual(name=legend_title, values=colors)
    } else {
        p1 <- p1 + scale_fill_discrete(name=legend_title)
    }
    
    # Add additional theme elements
    p1 <- p1 + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p1) }
    
    invisible(p1)
}