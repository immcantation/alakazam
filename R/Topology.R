# Ig lineage topology analysis

#### Graph analysis functions ####

#' Generate subtree summary statistics for a tree
#'
#' \code{summarizeSubtrees} calculates summary statistics for each node of a tree. Includes
#' both node properties and subtree properties.
#'
#' @param    graph   igraph object containing an annotated lineage tree.
#' @param    root    name of the root (germline) node.
#' @param    fields  vector of vertex annotation names to retain in the return data.frame.
#' 
#' @return   A data.frame with columns: 
#'           \itemize{
#'             \item  \code{NAME}:                node name
#'             \item  \code{PARENT}:              name of the parent node
#'             \item  \code{OUTDEGREE}:           number of edges leading from the node
#'             \item  \code{SUBTREE_SIZE}:        total number of nodes within the subtree rooted at the node
#'             \item  \code{SUBTREE_DEPTH}:       the depth of the subtree that is rooted at the node
#'             \item  \code{SUBTREE_PATHLENGTH}:  the maximum pathlength beneath the node
#'           }
#' 
#' @seealso  See \link{buildPhylipLineage} for generating input trees. 
#'           See \link{getPathLengths} for calculating path length to nodes.
#'           See \link{plotSubtrees} for...
#' 
#' @examples
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' V(graph)$label <- V(graph)$name
#'
#' # Plot graph with a tree layout
#' ly <- layout_as_tree(graph, root="Germline")
#' plot(graph, layout=ly)
#' 
#' # Summarize
#' summarizeSubtrees(graph, root="Germline")
#' 
#' @export
summarizeSubtrees <- function(graph, root="Germline", fields=NULL) {
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
    node_df$SUBTREE_SIZE <- apply(paths_step, 1, function(x) length(na.omit(x)))
    # Define number of levels below each node
    node_df$SUBTREE_DEPTH <- apply(paths_step, 1, max, na.rm=TRUE) + 1
    # Define the maximum shortest path length (genetic distance) to a leaf from each node
    node_df$SUBTREE_PATHLENGTH <- apply(paths_length, 1, max, na.rm=TRUE)
    
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
#'                     if \code{NULL} consider all nodes. This does not affect the weighted
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
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' V(graph)$label <- V(graph)$name
#'
#' # Plot graph with a tree layout
#' ly <- layout_as_tree(graph, root="Germline")
#' plot(graph, layout=ly)
#' 
#' # Consider all nodes
#' getPathLengths(graph, root="Germline")
#' 
#' # Exclude nodes without an isotype annotation from step count
#' getPathLengths(graph, root="Germline", field="isotype", exclude=NA)
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
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' V(graph)$label <- c("Germline", "Inferred", "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#'
#' # Plot graph with a tree layout
#' ly <- layout_as_tree(graph, root="Germline")
#' plot(graph, layout=ly)
#' 
#' # Use unweighted path length and do not exclude any nodes
#' getMRCA(graph, path="steps", root="Germline")
#'
#' # Exclude nodes without an isotype annotation and use weighted path length
#' getMRCA(graph, path="distance", root="Germline", field="isotype", exclude=NA)
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
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6, 6, 7, 7, 8, 7, 9))
#' V(graph)$name <- c("Germline", "Inferred1", "Seq1", "Seq2", "Seq3", 
#'                    "Seq4", "Inferred2", "Seq5", "Seq6")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA", NA, "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1, 1, 1, 1)
#' E(graph)$label <- E(graph)$weight
#' V(graph)$label <- V(graph)$isotype
#' 
#' # Plot graph with a tree layout
#' ly <- layout_as_tree(graph, root="Germline")
#' plot(graph, layout=ly)
#' 
#' # Count direct edges between isotypes including inferred nodes
#' tableEdges(graph, "isotype", exclude=NULL)
#' 
#' # Count direct edges excluding edges to and from inferred nodes
#' tableEdges(graph, "isotype", exclude=NA)
#' 
#' # Count indirect edges walking through inferred nodes
#' tableEdges(graph, "isotype", indirect=TRUE, exclude=NA)
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
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' V(graph)$label <- V(graph)$name
#' 
#' # Count edges between isotypes
#' permuteLabels(graph, "isotype")
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
#' @param    graphs   list of igraph object containing annotated lineage trees.
#' @param    field    string defining the annotation field to test.
#' @param    root     name of the root (germline) node.
#' @param    exclude  vector of strings defining \code{field} values to exclude from the
#'                    set of potential founder annotations.
#' @param    nperm    number of permutations to perform.
#' 
#' @return   A named list containing two data.frames summarizing the MRCA test (obs)
#'           and the null distribution of potential MRCAs (perm).
#'           
#' @seealso  Uses \link{getMRCA} and \link{getPathLengths}. Return object
#'           can be plotted with \link{plotMRCATest}.
#'           
#' @examples
#' # Define simple set of graphs
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' V(graph)$label <- V(graph)$name
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' graph2 <- graph
#' E(graph2)$weight <- c(10, 3, 3, 4, 1)
#' graph3 <- graph2
#' V(graph3)$isotype <- c(NA, NA, "IgM", "IgM", "IgA", "IgA")
#' 
#' # Count edges between isotypes
#' graphs <- list(A=graph, B=graph, C=graph2, D=graph3)
#' testMRCA(graphs, "isotype", nperm=100)
#' 
#' @export
testMRCA <- function(graphs, field, root="Germline", exclude=c("Germline", NA), 
                     nperm=200) {
    # TODO: should probably return a class instead of a list
    
    # Function to resolve ambiguous founders
    # @param  x      data.frame from getMRCA
    # @param  field  annotation field
    .resolveMRCA <- function(x, field) {
        x %>% filter_(interp(~!duplicated(x), x=as.name(field))) %>%
            filter_(interp(~length(x) == 1, x=as.name(field)))
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
            dplyr::summarize(COUNT=n()) %>%
            ungroup() %>%
            dplyr::mutate_(FREQ=interp(~x/sum(x, na.rm=TRUE), x=as.name("COUNT")))
        
        return(mrca_sum)
    }
    
    # Assign numeric names if graphs is an unnamed list
    if (is.null(names(graphs))) { names(graphs) <- 1:length(graphs) }
    
    # Summarize observed MRCA counts
    obs_sum <- .countMRCA(graphs, field=field, exclude=exclude)

    # Generate edge null distribution via permutation
    cat("-> PERMUTING TREES\n")
    pb <- txtProgressBar(min=0, max=nperm, initial=0, width=40, style=3)
    perm_list <- list()
    for (i in 1:nperm) {
        # Permute labels
        tmp_list <- lapply(graphs, permuteLabels, field=field, exclude=exclude)
        # Summarize MRCA counts
        tmp_sum <- .countMRCA(tmp_list, field=field, exclude=exclude)
        # Update permutation set
        tmp_sum$PERM <- i
        perm_list[[i]] <- tmp_sum
        
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    perm_sum <- bind_rows(perm_list)
    
    # Test observed against permutation distribution
    for (i in 1:nrow(obs_sum)) {
        x <- obs_sum$ANNOTATION[i]
        f <- ecdf(perm_sum$COUNT[perm_sum$ANNOTATION == x])
        obs_sum[i, "PVALUE"] <- 1 - f(obs_sum$COUNT[i])
    }
 
    return(list(obs=obs_sum, perm=perm_sum))
}


#' Tests for parent-child annotation enchrichment in lineage trees
#' 
#' \code{testEdges} performs a permutation test on a set of lineage trees to determine
#' the significance of an annotation's association with parent-child relationships.
#'
#' @param    graphs   list of igraph objects with vertex annotations.
#' @param    field    string defining the annotation field to permute.
#' @param    exclude  vector of strings defining \code{field} values to exclude from 
#'                    permutation.
#' @param    nperm    number of permutations to perform.
#' 
#' @return   A named list containing two data.frames summarizing the observed edge counts (obs)
#'           and the null distribution of edge counts (perm).
#' 
#' @seealso  Uses \link{tableEdges} and \link{permuteLabels}. Return object
#'           can be plotted with \link{plotEdgeTest}.
#'           
#' @examples
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' V(graph)$label <- V(graph)$name
#' 
#' # Count edges between isotypes
#' graphs <- list(A=graph, B=graph, C=graph)
#' testEdges(graphs, "isotype", nperm=100)
#' 
#' @export
testEdges <- function(graphs, field, exclude=c("Germline", NA), nperm=200) {
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
    cat("-> PERMUTING TREES\n")
    pb <- txtProgressBar(min=0, max=nperm, initial=0, width=40, style=3)
    perm_list <- list()
    for (i in 1:nperm) {
        # Permute annotations
        tmp_list <- lapply(graphs, permuteLabels, field=field, exclude=exclude)
        # Count edges
        tmp_sum <- .countEdges(tmp_list, field, exclude)
        # Update permutation set
        tmp_sum$PERM <- i
        perm_list[[i]] <- tmp_sum
        
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    perm_sum <- bind_rows(perm_list)
    
    # Test observed against permutation distribution
    for (i in 1:nrow(obs_sum)) {
        x <- obs_sum$PARENT[i]
        y <- obs_sum$CHILD[i]
        f <- ecdf(perm_sum$COUNT[perm_sum$PARENT == x & perm_sum$CHILD == y])
        obs_sum[i, "PVALUE"] <- 1 - f(obs_sum$COUNT[i])
    }
    
    return(list(obs=obs_sum, perm=perm_sum))
}


#### Plotting functions #####

#' Plot the results of an edge permutation test
#' 
#' \code{plotEdgeTest} plots the results of an edge permutation test performed with 
#' \code{testEdges}.
#'
#' @param    edge_test  named list returned by \code{testEdges}.
#' 
#' @return   NULL
#' 
#' @seealso  \link{testEdges}.
#' 
#' @examples
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' V(graph)$label <- V(graph)$name
#' 
#' # Count edges between isotypes
#' graphs <- list(A=graph, B=graph, C=graph)
#' edge_test <- testEdges(graphs, "isotype", nperm=100)
#' plotEdgeTest(edge_test)
#' 
#' @export
plotEdgeTest <- function(edge_test) {
    obs_sum <- edge_test[["obs"]]
    perm_sum <- edge_test[["perm"]]
    
    # Plot edge null distribution
    p1 <- ggplot(perm_sum, aes(x=count)) +
        ggtitle("Edge Null Distribution") +
        getBaseTheme() + 
        xlab("Number of edges") +
        ylab("Number of realizations") + 
        geom_histogram(fill="steelblue", color="black", binwidth=1) +
        #geom_density(adjust=3.0, color="dimgrey", size=1.0, linetype=1) + 
        geom_vline(data=obs_sum, aes(xintercept=count), color="firebrick", linetype=3, size=1.25) + 
        facet_grid(child ~ parent, labeller=label_both, scales="free")
    plot(p1)
    
    # Plot ECDF of edge null distribution
    p2 <- ggplot(perm_sum, aes(x=count)) +
        ggtitle("Edge ECDF") +
        getBaseTheme() + 
        xlab("Number of edges") +
        ylab("P-Value") +
        stat_ecdf(color="black", size=1) +
        geom_vline(data=obs_sum, aes(xintercept=count), color="firebrick", linetype=3, size=1.25) + 
        geom_hline(data=obs_sum, aes(yintercept=pvalue), color="steelblue", linetype=3, size=1.25) + 
        facet_grid(child ~ parent, labeller=label_both, scales="free")
    plot(p2)
}


#' Plot the results of a founder permutation test
#' 
#' \code{plotFounderTest} plots the results of a founder permutation test performed with 
#' \code{testFounders}.
#'
#' @param    founder_test  named list returned by \code{testFounders}.
#' 
#' @return   NULL
#' 
#' @seealso  \link{testFounders}.
#' 
#' @examples
#' # Define simple graph
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' V(graph)$label <- V(graph)$name
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' 
#' # Count edges between isotypes
#' graph2 <- graph
#' E(graph2)$weight <- c(10, 3, 3, 4, 1)
#' graph3 <- graph2
#' V(graph3)$isotype <- c(NA, NA, "IgM", "IgM", "IgA", "IgA")
#' graphs <- list(A=graph, B=graph2, C=graph3)
#' founder_test <- testFounders(graphs, "isotype", nperm=100)
#' 
#' # Plot
#' plotFounderTest(founder_test)
#' 
#' @export
plotFounderTest <- function(founder_test) {
    obs_sum <- founder_test[["obs"]]
    perm_sum <- founder_test[["perm"]]
    
    # Plot founder null distribution
    p1 <- ggplot(perm_sum, aes(x=count)) +
        ggtitle("Founder Null Distribution") +
        getBaseTheme() + 
        xlab("Number of founders") +
        ylab("Number of realizations") + 
        geom_histogram(fill="steelblue", color="black", binwidth=1) +
        #geom_density(adjust=3.0, color="dimgrey", size=1.0, linetype=1) + 
        geom_vline(data=obs_sum, aes(xintercept=count), color="firebrick", linetype=3, size=1.25) + 
        facet_wrap(~ label, ncol=1, scales="free_y")
    plot(p1)
    
    # Plot ECDF of founder null distribution
    p2 <- ggplot(perm_sum, aes(x=count)) +
        ggtitle("Founder ECDF") +
        getBaseTheme() + 
        xlab("Number of founders") +
        ylab("P-Value") +
        stat_ecdf(color="black", size=1) +
        geom_vline(data=obs_sum, aes(xintercept=count), color="firebrick", linetype=3, size=1.25) + 
        geom_hline(data=obs_sum, aes(yintercept=pvalue), color="steelblue", linetype=3, size=1.25) + 
        facet_wrap(~ label, nrow=1, scales="free_y")
    plot(p2)
}


#' Plots subtree properties distributions for multiple trees
#' 
#' \code{plotsSubtrees} summarizes the subtree sizes, depths and pathlengths by vertex
#' annotation values and plots them.
#'
#' @param    graphs   list of igraph objects with vertex annotations.
#' @param    field    string defining the annotation field.
#' @param    root     name of the root (germline) node.
#' @param    exclude  vector of strings defining \code{field} values to exclude from the results.
#' @param    style    string specifying the style of plot to draw. One of c("box", "violin")
#' 
#' @return   A data.frame of the subtree properties for all trees.
#' 
#' @seealso  \link{summarizeSubtrees}.
#' 
#' @examples
#' # Define simple set of graphs
#' library(igraph)
#' graph <- make_directed_graph(c(1, 2, 2, 3, 2, 4, 3, 5, 3, 6))
#' V(graph)$name <- c("Germline", "Inferred", "Seq1", "Seq2", "Seq3", "Seq4")
#' V(graph)$isotype <- c(NA, NA, "IgM", "IgG", "IgA", "IgA")
#' V(graph)$label <- V(graph)$name
#' E(graph)$weight <- c(10, 3, 6, 4, 1)
#' E(graph)$label <- E(graph)$weight
#' graph2 <- graph
#' E(graph2)$weight <- c(10, 3, 3, 4, 1)
#' graph3 <- graph2
#' V(graph3)$isotype <- c(NA, NA, "IgM", "IgM", "IgA", "IgA")
#' 
#' # Plot subtrees
#' graphs <- list(A=graph, B=graph, C=graph2, D=graph3)
#' df <- plotSubtrees(graphs, "isotype")
#' 
#' @export
plotSubtrees <- function(graphs, field, root="Germline", exclude=c("Germline", NA), style="box") {
    if (!(style %in% c("box", "violin"))) {
        stop("The style argument must be one of c('box', 'violin')")
    }
    
    # Assign numeric names if graphs is an unnamed list
    if (is.null(names(graphs))) { names(graphs) <- 1:length(graphs) }
    
    # Get subtree summarizes
    sum_df <- plyr::ldply(graphs, summarizeSubtrees, root=root, fields=field, .id="tree_identifier")
    sum_df <- sum_df[!(sum_df[, field] %in% exclude), ] 
    
    # Normalize subtree properties within tree
    sum_df <- plyr::ddply(sum_df, c("tree_identifier"), plyr::mutate, 
                          outdegree_norm=outdegree/length(name),
                          size_norm=subtree_size/max(subtree_size, na.rm=TRUE),
                          depth_norm=subtree_depth/max(subtree_depth, na.rm=TRUE),
                          pathlength_norm=subtree_pathlength/max(subtree_pathlength, na.rm=TRUE))
        
    # Plot outdegree
    p1 <- ggplot(subset(sum_df, is.finite(outdegree_norm)), aes_string(x=field, y="outdegree_norm")) + 
        ggtitle("Outdegree") + 
        getBaseTheme() + theme(legend.position="none") +
        xlab("") +
        ylab("Outdegree (percent of tree size)") +
        scale_y_continuous(labels=percent) +
        expand_limits(y=0)
    if (style == "box") {
        p1 <- p1 + geom_boxplot(aes_string(fill=field), width=0.7, alpha=0.8)
    } else if (style == "violin") {
        p1 <- p1 + geom_violin(aes_string(fill=field), adjust=1.5, scale="width", trim=T, width=0.7, alpha=0.8) +
            geom_errorbar(stat="hline", yintercept="mean", aes(ymin=..y.., ymax=..y..),
                          width=0.8, size=2.5, color="black")
    }
    plot(p1)

    # Plot subtree size
    p2 <- ggplot(subset(sum_df, is.finite(size_norm)), aes_string(x=field, y="size_norm")) + 
        ggtitle("Subtree size") + 
        getBaseTheme() + theme(legend.position="none") +
        xlab("") +
        ylab("Size (percent of tree size)") +
        scale_y_continuous(labels=percent) +
        expand_limits(y=0)
    if (style == "box") {
        p2 <- p2 + geom_boxplot(aes_string(fill=field), width=0.7, alpha=0.8)
    } else if (style == "violin") {
        p2 <- p2 + geom_violin(aes_string(fill=field), adjust=1.5, scale="width", trim=T, width=0.7, alpha=0.8) +
            geom_errorbar(stat="hline", yintercept="mean", aes(ymin=..y.., ymax=..y..),
                          width=0.8, size=2.5, color="black")
    }
    plot(p2)
    
    # Plot normalized depth of subtree under node
    p3 <- ggplot(subset(sum_df, is.finite(depth_norm)), aes_string(x=field, y="depth_norm")) + 
        ggtitle("Subtree depth") + 
        getBaseTheme() + theme(legend.position="none") +
        xlab("") +
        ylab("Depth (percent of tree depth)") +
        scale_y_continuous(labels=percent) +
        expand_limits(y=0)
    if (style == "box") {
        p3 <- p3 + geom_boxplot(aes_string(fill=field), width=0.7, alpha=0.8)
    } else if (style == "violin") {
        p3 <- p3 + geom_violin(aes_string(fill=field), adjust=1.5, scale="width", trim=T, width=0.7, alpha=0.8) +
            geom_errorbar(stat="hline", yintercept="mean", aes(ymin=..y.., ymax=..y..),
                          width=0.8, size=2.5, color="black")
    }
    plot(p3)
 
    # Plot normalized pathlength of subtree under node
    p4 <- ggplot(subset(sum_df, is.finite(pathlength_norm)), aes_string(x=field, y="pathlength_norm")) + 
        ggtitle("Subtree maximum pathlength to leaf") + 
        getBaseTheme() + theme(legend.position="none") +
        xlab("") +
        ylab("Maximum pathlength (percent of tree pathlength)") +
        scale_y_continuous(labels=percent) +
        expand_limits(y=0)
    if (style == "box") {
        p4 <- p4 + geom_boxplot(aes_string(fill=field), width=0.7, alpha=0.8)
    } else if (style == "violin") {
        p4 <- p4 + geom_violin(aes_string(fill=field), adjust=1.5, scale="width", trim=T, width=0.7, alpha=0.8) +
            geom_errorbar(stat="hline", yintercept="mean", aes(ymin=..y.., ymax=..y..),
                          width=0.8, size=2.5, color="black")
    }
    plot(p4)
    
    return(sum_df)
}