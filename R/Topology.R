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
#'           See \link{getPathLengths} and \link{plotSubtrees} for...
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
#'           \link{getFounder}.
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
    # TODO:  should steps be level?
    # Define path length data.frame
    path_df <- data.frame(name=V(graph)$name, stringsAsFactors=FALSE)

    # Get indices of excluded vertices
    skip_idx <- which(path_df$name == root)
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
        path_df[i, "steps"] <- sum(!(v %in% skip_idx)) 
        path_df[i, "distance"] <- sum(E(graph, path=v)$weight)
    }
    
    return(path_df)
}


#' Retrieve the first non-root node of a lineage tree
#' 
#' \code{getFounder} returns the set of founding node(s) for an igraph object. The founding 
#' node(s) are defined as the set of nodes within the minimum weighted or unweighted path length
#' from the root (germline) of the lineage tree, allowing for exclusion of specific groups of
#' nodes.
#'
#' @param    graph    igraph object with vertex annotations.
#' @param    path     string defining whether to use unweighted (steps) or weighted (distance) 
#'                    measures for determining the founder node set. Must be one of 
#'                    \code{c("steps", "distance")}. 
#' @param    root     name of the root (germline) node.
#' @param    field    annotation field to use for both unweighted path length exclusion and
#'                    consideration as a founder node. if \code{NULL} do not exclude any nodes.
#' @param    exclude  vector of annotation values in the given field to exclude from potential 
#'                    founder set. If \code{NULL} do not exclude any nodes. Has no effect if 
#'                    \code{field=NULL}.
#'                    
#' @return   A data.frame of the founder node(s).
#' 
#' @seealso  \link{getPathLengths}.
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
#' getFounder(graph, path="steps", root="Germline")
#'
#' # Exclude nodes without an isotype annotation and use weighted path length
#' getFounder(graph, path="distance", root="Germline", field="isotype", exclude=NA)
#' 
#' @export
getFounder <- function(graph, path="distance", root="Germline", field=NULL, exclude=NULL) {
    # Get distance from root
    path_df <- getPathLengths(graph, root=root, field=field, exclude=exclude)
    
    # Get indices of excluded vertices
    skip_idx <- which(path_df$name == root)
    if (!is.null(field)) {
        g <- vertex_attr(graph, name=field)
        skip_idx <- union(skip_idx, which(g %in% exclude))
    }
    
    # Get founder nodes
    if (path == "distance") { 
        path_len <- setNames(path_df$distance, 1:nrow(path_df))
    } else if (path == "steps") {
        path_len <- setNames(path_df$steps, 1:nrow(path_df))
    } else {
        stop("Invalid value for 'path' parameter. Must be one of c('distance', 'steps').\n")
    }
    path_len <- path_len[-skip_idx]
    root_idx <- as.numeric(names(path_len)[which(path_len == min(path_len))])
    root_df <- igraph::as_data_frame(graph, what="vertices")[root_idx, ]
    root_df$steps <- path_df$steps[root_idx]
    root_df$distance <- path_df$distance[root_idx]
    
    return(root_df)
}


#' Tabulate the number of edges between annotations within a tree
#'
#' \code{tableEdges} creates a table of the total number of connections (edges) for each 
#' unique pair of annotations within the tree over all vertices.
#' 
#' @param    graph     igraph object with vertex annotations.
#' @param    field     string defining the annotation field to count.
#' @param    indirect  if \code{FALSE} count direct connections (edges) only. If 
#'                     \code{TRUE} walk through any nodes with annotations specified in 
#'                     the \code{argument} to count indirect connections. Specifying
#'                     \code{indirect=TRUE} with \code{exclude=NULL} will have no effect.
#' @param    exclude   vector of strings defining \code{field} values to exclude from counts.
#'                     Edges that either start or end with the specified annotations will not
#'                     be counted. If \code{NULL} count all edges.
#'                     
#' @return   A data.frame with columns (parent, child, count) defining total
#'           annotation connections in the tree.
#' 
#' @seealso  link{\code{testEdges}}.
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
                edge_list[[i]] <- data.frame("parent"=parent, "child"=children)
            }
        }
        
        # Merge edge list into data.frame
        edge_df <- plyr::rbind.fill(edge_list)        
    }
    else {
        # Get adjacency list
        edge_mat <- as_edgelist(graph, names=FALSE)
        edge_mat <- vertex_attr(graph, name=field, index=edge_mat)
        edge_mat <- matrix(edge_mat, ncol=2, dimnames=list(NULL, c("parent", "child")))

        # Build and subset edge data.frame
        edge_df <- as.data.frame(edge_mat, stringsAsFactors=FALSE)
        edge_df <- subset(edge_df, !(parent %in% exclude) & !(child %in% exclude))
    }
    
    # Count edges
    edge_tab <- plyr::ddply(edge_df, c("parent", "child"), summarize, count=length(parent))
    
    return(edge_tab)
}


#' Permute the node labels of a tree
#' 
#' \code{permuteLabels} permutes the node annotations of a lineage tree.
#'
#' @param    graph    igraph object with vertex annotations.
#' @param    field    string defining the annotation field to permute.
#' @param    exclude  vector of strings defining \code{field} values to exclude from permutation.
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


#' Infers the founding node annotation of a list of trees
#'
#' @param    graphs   list of igraph objects with vertex annotations.
#' @param    field    string defining the annotation field to test.
#' @param    root     name of the root (germline) node.
#' @param    exclude  vector of strings defining \code{field} values to exclude from the
#'                    set of potential founder annotations.
#' @param    nperm    number of permutations to perform.
#' 
#' @return   A named list containing two data.frames summarizing the founder test (obs)
#'           and the null distribution of potential founders (perm).
#'           
#' @seealso  Uses \link{getFounder} and \link{getPathLengths}. Return object
#'           can be plotted with \link{plotFounderTest}.
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
#' testFounders(graphs, "isotype", nperm=100)
#' 
#' @export
testFounders <- function(graphs, field, root="Germline", exclude=c("Germline", NA), 
                         nperm=2000) {
    # Function to resolve ambiguous founders
    resolveFounders <- function(df) {
        # Remove duplicate labels
        df <- plyr::ddply(df, c("graph"), subset, !duplicated(label))
        # Discard unresolved ties
        df <- plyr::ddply(df, c("graph"), subset, length(label) == 1)
        
        return(df)
    }
    
    # Assign numeric names if graphs is an unnamed list
    if (is.null(names(graphs))) { names(graphs) <- 1:length(graphs) }
        
    # Get founders
    obs_df <- plyr::ldply(graphs, getFounder, path="distance", field=field, 
                    exclude=exclude, .id="graph")
    obs_df <- obs_df[, c("graph", field)]
    names(obs_df) <- c("graph", "label")
    obs_df <- resolveFounders(obs_df)
    # Summarize founder counts  
    obs_sum <- plyr::ddply(obs_df, c("label"), plyr::summarize, count=length(label))
    obs_sum <- transform(obs_sum, freq=count/sum(count, na.rm=TRUE))
    
    # Generate edge null distribution via permutation
    cat("-> PERMUTING TREES\n")
    pb <- txtProgressBar(min=0, max=nperm, initial=0, width=40, style=3)
    perm_list <- list()
    for (i in 1:nperm) {
        tmp_list <- plyr::llply(graphs, permuteLabels, field=field, exclude=exclude)
        tmp_df <- plyr::ldply(tmp_list, getFounder, path="distance", field=field, 
                        exclude=exclude, .id="graph")
        tmp_df <- tmp_df[, c("graph", field)]
        names(tmp_df) <- c("graph", "label")
        tmp_df <- resolveFounders(tmp_df)
        tmp_sum <- plyr::ddply(tmp_df, c("label"), plyr::summarize, count=length(label))
        tmp_sum$perm <- i
        perm_list[[i]] <- tmp_sum
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    perm_sum <- plyr::rbind.fill(perm_list)
    
    # Test observed against permutation distribution
    for (i in 1:nrow(obs_sum)) {
        x <- obs_sum$label[i]
        f <- ecdf(perm_sum$count[perm_sum$label == x])
        obs_sum[i, "pvalue"] <- 1 - f(obs_sum$count[i])
    }
 
    return(list(obs=obs_sum, perm=perm_sum))
}


#' Tests for edge enrichment via permutation
#' 
#' \code{testEdges} tests the significance of edge abundance between vertext annotation pairs 
#' within a list of trees via permutation of the graph labels.
#'
#' @param    graphs   list of igraph objects with vertex annotations.
#' @param    field    string defining the annotation field to permute.
#' @param    exclude  vector of strings defining \code{field} values to exclude from permutation.
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
testEdges <- function(graphs, field, exclude=c("Germline", NA), nperm=2000) {
    # Assign numeric names if graphs is an unnamed list
    #if (is.null(names(graphs))) { names(graphs) <- 1:length(graphs) }
    
    # Count edges of unpermuted data
    obs_df <- plyr::ldply(graphs, tableEdges, field=field, exclude=exclude)
    obs_sum <- plyr::ddply(obs_df, c("parent", "child"), plyr::summarize, count=sum(count, na.rm=TRUE))

    # Generate edge null distribution via permutation
    cat("-> PERMUTING TREES\n")
    pb <- txtProgressBar(min=0, max=nperm, initial=0, width=40, style=3)
    perm_list <- list()
    for (i in 1:nperm) {
        tmp_list <- plyr::llply(graphs, permuteLabels, field=field, exclude=exclude)
        tmp_df <- plyr::ldply(tmp_list, tableEdges, field=field, exclude=exclude)
        tmp_sum <- plyr::ddply(tmp_df, c("parent", "child"), 
                               plyr::summarize, count=sum(count, na.rm=TRUE))
        tmp_sum$perm <- i
        perm_list[[i]] <- tmp_sum
        setTxtProgressBar(pb, i)
    }
    cat("\n")
    perm_sum <- plyr::rbind.fill(perm_list)
    
    # Test observed against permutation distribution
    for (i in 1:nrow(obs_sum)) {
        x <- obs_sum$parent[i]
        y <- obs_sum$child[i]
        f <- ecdf(perm_sum$count[perm_sum$parent == x & perm_sum$child == y])
        obs_sum[i, "pvalue"] <- 1 - f(obs_sum$count[i])
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