# Ig lineage reconstruction via maximum parsimony
# 
# @author     Jason Anthony Vander Heiden
# @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2014.12.10


#### Classes ####

#' S4 class defining a clone
#' 
#' \code{ChangeoClone} defines common fields to perform lineage recontruction
#' from a Change-O clone.
#' 
#' @slot  data        data.frame containing SEQUENCE_ID, SEQUENCE_GAP,
#'                    and additional annotation fields.
#' @slot  clone       string defining the clone identifier.
#' @slot  germline    string containing the germline sequence for the clone.
#' @slot  v_gene      string defining the V segment gene call.
#' @slot  v_gene      string defining the J segment gene call.
#' @slot  junc_len    numeric junction length.
#' 
#' @name ChangeoClone
#' @export
setClass("ChangeoClone", 
         slots=c(data="data.frame",
                 clone="character",
                 germline="character", 
                 v_gene="character", 
                 j_gene="character", 
                 junc_len="numeric"))


#### Preprocessing functions ####

#' Process a clone for lineage construction
#' 
#' \code{prepChangeoClone} preprocessing a clonal group for lineage reconstruction.
#' 
#' Takes a data.frame with Change-O style column as input and masks gap positions, 
#' masks ragged ends, remove duplicates sequences, and merges annotations associated 
#' with duplicate sequences.
#' 
#' The clone identifier, germline sequence, V gene, J gene, and junction length are 
#' determined from the first entry in the CLONE, GERMLINE_GAP_D_MASK, V_CALL, J_CALL, 
#' and JUNCTION_GAP_LENGTH columns, respectively. For any given clone, each value in these
#' columns should be identical.
#'
#' @param    data         a data.frame containing the Change-O data for a clone. 
#'                        The data.frame must contain the following columns:
#'                        SEQUENCE_ID, SEQUENCE_GAP, CLONE, GERMLINE_GAP_D_MASK
#'                        V_CALL, J_CALL, JUNCTION_GAP_LENGTH.
#' @param    max_mask     maximum number of characters to mask at the leading and trailing
#'                        sequence ends. If \code{NULL} then the threshold is automatically
#'                        set to the highest number of observed outer ends.  If set to 0 then
#'                        no masking is performed.
#' @param    text_fields  text annotation columns to process during duplicate removal.
#' @param    num_fields   numeric annotation columns to process during duplicate removal.
#' @return   a \code{ChangeoClone} object containing the modified clone
#' 
#' @seealso  Executes in order \code{\link{maskSeqGaps}}, \code{\link{maskSeqEnds}}
#'           and \code{\link{collapseDuplicates}}. Returns a \code{\link{ChangeoClone}}.
#' @examples
#' # Example Change-O data.frame
#' df <- data.frame(SEQUENCE_ID=LETTERS[1:4],
#'                  SEQUENCE_GAP=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
#'                  V_CALL="Homsap IGKV1-39*01 F",
#'                  J_CALL="Homsap IGKJ5*01 F",
#'                  JUNCTION_GAP_LENGTH=2,
#'                  GERMLINE_GAP_D_MASK="CCCCAGGG",
#'                  CLONE=1,
#'                  TYPE=c("IgM", "IgG", "IgG", "IgA"),
#'                  COUNT=1:4,
#'                  stringsAsFactors=FALSE)
#' 
#' # Without end masking
#' prepChangeoClone(df, text_fields="TYPE", num_fields="COUNT")
#'
#' # With end masking
#' prepChangeoClone(df, max_mask=3, text_fields="TYPE", num_fields="COUNT")
#'
#' @export
prepChangeoClone <- function(data, max_mask=0, text_fields=NULL, num_fields=NULL) {
    # Replace gaps with Ns and masked ragged ends
    tmp_df <- data[, c("SEQUENCE_ID", "SEQUENCE_GAP", text_fields, num_fields)]
    tmp_df[, "SEQUENCE_GAP"] <- maskSeqGaps(tmp_df[, "SEQUENCE_GAP"], outer_only=FALSE)
    tmp_df[, "SEQUENCE_GAP"] <- maskSeqEnds(tmp_df[, "SEQUENCE_GAP"], max_mask=max_mask, trim=FALSE)
    
    # Remove duplicates
    tmp_df <- collapseDuplicates(tmp_df, id="SEQUENCE_ID", seq="SEQUENCE_GAP",
                                 text_fields=text_fields, num_fields=num_fields)
    
    # Define return object
    clone <- new("ChangeoClone", 
                 data=tmp_df,
                 clone=as.character(data[1, "CLONE"]),
                 germline=maskSeqGaps(data[1, "GERMLINE_GAP_D_MASK"], outer_only=FALSE), 
                 v_gene=getGene(data[1, "V_CALL"]), 
                 j_gene=getGene(data[1, "J_CALL"]), 
                 junc_len=data[1, "JUNCTION_GAP_LENGTH"])
    
    return(clone)
}


#### PHYLIP functions ####

# Create PHYLIP input files in a temporary folder
#
# @param   clone  a ChangeoClone object
# @param   path   a directory to store the write the output files to
# @return  a named vector translating SEQUENCE_ID (names) to PHYLIP taxa (values)
writePhylipInput <- function(clone, path) {
    # Define PHYLIP columns
    nseq <- nrow(clone@data)
    v1 <- c(sprintf('%-9s', nseq + 1),
            sprintf("%-9s", "Germline"), 
            sprintf("SAM%-6s", 1:nseq))
    v2 <- c(nchar(clone@germline),
            clone@germline, 
            clone@data[, "SEQUENCE_GAP"])
    phy_df <- data.frame(v1, v2, stringsAsFactors=F)
    
    # Define names vector mapping taxa names to original sequence identifiers
    id_map <- setNames(str_trim(v1[-(1:2)]), clone@data[, "SEQUENCE_ID"])
    
    # Create PHYLIP input file
    write.table(phy_df, file=file.path(path, "infile"), 
                quote=F, sep=" ", col.names=F, row.names=F)    
    
    return(id_map)
}


# Run PHYLIP dnapars application
#
# @param   path          temporary directory containing infile.
# @param   dnapars_exec  path to the dnapars executable.
# @param   verbose       if TRUE suppress phylip console output.
# @return  NULL
runPhylip <- function(path, dnapars_exec, verbose=FALSE) {
    # Remove old files
    if (file.exists(file.path(path, "outfile"))) { file.remove(file.path(path, "outfile")) }
    if (file.exists(file.path(path, "outtree"))) { file.remove(file.path(path, "outtree")) }    

    # Run dnapars
    phy_options <- c("S", "Y", "I", "4", "5", ".")
    if (verbose) {
        system(paste0("cd ", path, "; ", dnapars_exec), input=c(phy_options, "Y"))
    } else {
        system(paste0("cd ", path, "; ", dnapars_exec), input=c(phy_options, "Y"),
               ignore.stdout=TRUE, ignore.stderr=TRUE)
    }
}


# Reads in the PHYLIP outfile
#
# @param   path  the temporary folder containing the dnapars outfile
# @return  a character vector with each item as a line in the outfile
readPhylipOutput <- function(path) {
    phylip_out <- scan(file.path(path, "outfile"), what="character", sep="\n", 
                       blank.lines.skip=FALSE, strip.white=FALSE, quiet=TRUE)
    return(phylip_out)
}


# Test for successful PHYLIP dnapars run by checking the outfile
#
# @param   phylip_out  a character vector returned by readPhylipOut
# @return  TRUE if trees built 
#          FALSE if no trees built
checkPhylipOutput <- function(phylip_out) {
    # Check for failed tree build
    result <- !(any(grepl('-1 trees in all found', phylip_out)))
    
    return(result)
}


# Extracts inferred sequences from PHYLIP dnapars outfile
#
# @param   phylip_out   a character vector returned by readPhylipOutput
# @return  a list containing an id vector, a sequence vector and an annotation data.frame
getPhylipInferred <- function(phylip_out) {
    # Process dnapars output
    seq_start <- min(grep("From\\s+To\\s+Any Steps\\?\\s+State at upper node", 
                          phylip_out, perl=T, fixed=F))
    seq_empty <- grep("^\\s*$", phylip_out[seq_start:length(phylip_out)], perl=T, fixed=F)
    seq_len <- seq_empty[min(which(seq_empty[-1] == (seq_empty[-length(seq_empty)] + 1)))]
    seq_block <- paste(phylip_out[(seq_start + 2):(seq_start + seq_len - 2)], collapse="\n")
    seq_df <- read.table(textConnection(seq_block), as.is=T, fill=T, blank.lines.skip=F)
    
    # Correct first line of block and remove blank rows
    fix.row <- c(1, which(is.na(seq_df[,1])) + 1)
    end_col <-  ncol(seq_df) - 2
    #seq_df[fix.row, ] <- cbind(0, seq_df[fix.row, 1], "no", seq_df[fix.row, 2:5], stringsAsFactors=F)
    seq_df[fix.row, ] <- cbind(0, seq_df[fix.row, 1], "no", seq_df[fix.row, 2:end_col], stringsAsFactors=F)
    seq_df <- seq_df[-(fix.row[-1] - 1), ]
    
    # Create data.frame of inferred sequences
    inferred_num <- unique(grep("^[0-9]+$", seq_df[, 2], value=T))
    inferred_seq <- sapply(inferred_num, function(n) { paste(t(as.matrix(seq_df[seq_df[, 2] == n, -c(1:3)])), collapse="") })
    
    return(data.frame(SEQUENCE_ID=paste0("Inferred", inferred_num), SEQUENCE_GAP=inferred_seq))
}


# Extracts graph edge list from a PHYLIP dnapars outfile
#
# @param   phylip_out  character vector returned by readPhylipOutput
# @param   id_map      named vector of PHYLIP taxa names (values) to sequence 
#                      identifiers (names) that will be translated. If NULL
#                      no taxa name translation is performed
# @return  a data.frame of edges with columns (from, to, weight)
getPhylipEdges <- function(phylip_out, id_map=NULL) {
    # Process dnapars output
    edge_start <- min(grep('between\\s+and\\s+length', phylip_out, 
                           perl=TRUE, fixed=FALSE))
    edge_len <- min(grep('^\\s*$', phylip_out[edge_start:length(phylip_out)], 
                         perl=TRUE, fixed=FALSE))
    edge_block <- paste(phylip_out[(edge_start + 2):(edge_start + edge_len - 2)], collapse='\n')
    edge_df <- read.table(textConnection(edge_block), col.names=c('from', 'to', 'weight'), 
                          as.is=TRUE)

    # Modify inferred taxa names to include "Inferred"
    inf_map <- unique(grep("^[0-9]+$", c(edge_df$from, edge_df$to), value=T))
    names(inf_map) <- paste0("Inferred", inf_map)
    edge_df$from <- translateStrings(edge_df$from, inf_map)
    edge_df$to <- translateStrings(edge_df$to, inf_map)
    
    if (!is.null(id_map)) {
        # Reassign PHYLIP taxa names to sequence IDs
        edge_df$from <- translateStrings(edge_df$from, id_map)
        edge_df$to <- translateStrings(edge_df$to, id_map)
    }
    
    return(edge_df)
}


# Modify edges of phylip output
#
# @param   edges    data.frame of edges returned by getPhylipEdges
# @param   clone    a ChangeoClone object containg sequence data
# @param   nuc_mat  nucleotide character distance matrix
# @return  a list of modified edges data.frame and clone object
modifyPhylipEdges <- function(edges, clone, nuc_mat=getNucMatrix(gap=0)) {
    # Move germline to root position
    germ_idx <- which(edges$to == "Germline")
    edges[germ_idx, c('from', 'to')] <- edges[germ_idx, c('to', 'from')]
    
    # Calculate edge mutations
    for (i in 1:nrow(edges)) {
        if (edges$from[i] == "Germline") {
            seq1 <- clone@germline
        } else {
            seq1 <- clone@data[clone@data[, "SEQUENCE_ID"] == edges$from[i], "SEQUENCE_GAP"]
        }
        seq2 <- clone@data[clone@data[, "SEQUENCE_ID"] == edges$to[i], "SEQUENCE_GAP"]
        edges$weight[i] <- getSeqDistance(seq1, seq2, nuc_mat)        
    }
    
    # Find rows zero weight edges with inferred parent nodes
    remove_row <- which(edges$weight == 0 & 
                        edges$from != "Germline" & 
                        grepl('^Inferred\\d+$', edges$from))
    
    # Replace inferred parent nodes with child nodes when edge weight is zero
    while (length(remove_row) > 0) {
        # Remove first node with zero distance to parent
        r <- remove_row[1]
        r_idx <- which(edges[c('from', 'to')] == edges$from[r], arr.ind=T)
        edges[r_idx] <- edges$to[r]
        
        # Recalculate edge weights for modified rows
        r_mod <- r_idx[, 1][r_idx[, 1] != r]
        for (i in r_mod) {
            if (edges$from[i] == "Germline") {
                seq1 <- clone@germline
            } else {
                seq1 <- clone@data[clone@data[, "SEQUENCE_ID"] == edges$from[i], "SEQUENCE_GAP"]
            }
            seq2 <- clone@data[clone@data[, "SEQUENCE_ID"] == edges$to[i], "SEQUENCE_GAP"]
            edges$weight[i] <- getSeqDistance(seq1, seq2, nuc_mat)      
        }
        
        # Remove row
        edges <- edges[-r, ]
        
        # Re-determine rows to remove
        remove_row <- which(edges$weight == 0 & 
                            edges$from != "Germline" & 
                            grepl('^Inferred\\d+$', edges$from))      
    }
    
    # Remove rows from clone
    keep_clone <- clone@data[, "SEQUENCE_ID"] %in% unique(c(edges$from, edges$to))
    clone@data <- clone@data[keep_clone, ]
    
    return(list(edges=edges, clone=clone))
}

# Convert edge data.frame and clone object to igraph graph object
#
# @param   edges  data.frame of edges returned by getPhylipEdges
# @param   clone  a ChangeoClone object containg sequence data
# @return  an igraph graph object
phylipToGraph <- function(edges, clone) {
    # Create igraph object
    g <- graph.data.frame(edges, directed=T)
    
    # Add germline sequence
    germ_idx <- which(V(g)$name == "Germline")
    g <- set.vertex.attribute(g, "sequence", index=germ_idx, clone@germline)
    
    # Add sample sequences and names
    clone_idx <- match(clone@data[, "SEQUENCE_ID"], V(g)$name) 
    g <- set.vertex.attribute(g, "sequence", index=clone_idx, clone@data[, "SEQUENCE_GAP"])
    
    # Add annotations
    ann_fields <- names(clone@data)[!(names(clone@data) %in% c("SEQUENCE_ID", "SEQUENCE_GAP"))]
    for (n in ann_fields) {
        g <- set.vertex.attribute(g, n, index=germ_idx, NA)
        g <- set.vertex.attribute(g, n, index=clone_idx, clone@data[, n])
    }
    
    # Add edge and vertex labels
    V(g)$label <- V(g)$name
    E(g)$label <- E(g)$weight
    
    # Add graph attributes
    g$clone <- clone@clone
    g$v_gene <- clone@v_gene
    g$j_gene <- clone@j_gene
    g$junc_len <- clone@junc_len
    
    return(g)
}


#' Infer an Ig lineage using PHYLIP
#' 
#' \code{buildPhylipLineage} reconstructs an Ig lineage using the dnapars
#' application of the PHYLIP package.
#'
#' @param    clone         \code{ChangeoClone} object containg clone data.
#' @param    dnapars_exec  path to the PHYLIP dnapars executable.
#' @param    rm_temp       if TRUE delete the temporary directory after running PHYLIP;
#'                         if FALSE keep the temporary directory.
#' @param    verbose       if \code{FALSE} suppress the output of dnapars; 
#'                         if \code{TRUE} STDOUT and STDERR of dnapars will be passed to the console.                        
#' @return   an igraph \code{graph} object
#' 
#' @seealso See \code{\link{igraph}} and \code{\link{igraph.plotting}} for working 
#'          with igraph \code{graph} objects.
#' @examples
#' \dontrun{
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Preprocess clone
#' clone <- subset(df, CLONE == 164)
#' clone <- prepChangeoClone(clone, text_fields=c("SAMPLE", "ISOTYPE"), num_fields="DUPCOUNT")
#' 
#' # Run PHYLIP and process output
#' dnapars_exec <- "~/apps/phylip-3.69/dnapars"
#' graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)
#' 
#' # Plot graph with a tree layout
#' ly <- layout.reingold.tilford(graph, root="Germline", circular=F, flip.y=T)
#' plot(graph, layout=ly)
#' }
#' 
#' @export
buildPhylipLineage <- function(clone, dnapars_exec, rm_temp=FALSE, verbose=FALSE) {
    if (nrow(clone@data) < 2) {
        warning("Clone ", clone@clone, " was skipped as it does not contain at least 
                2 unique sequences")
        return(NULL)
    }
    # Create temporary directory
    temp_path <- makeTempDir(paste0(clone@clone, "-phylip"))
    
    # Run PHYLIP
    id_map <- writePhylipInput(clone, temp_path)
    runPhylip(temp_path, dnapars_exec, verbose=verbose)
    phylip_out <- readPhylipOutput(temp_path)
    
    # Remove temporary directory
    if (rm_temp) {
        unlink(temp_path, recursive=TRUE)
    }
    
    # Check output for trees
    if (!checkPhylipOutput(phylip_out)) {
        warning('PHYLIP failed to generate trees for clone ', clone)
        return(NULL)
    }
    
    # Extract inferred sequences from PHYLIP output
    inf_df <- getPhylipInferred(phylip_out)
    clone@data <- rbind.fill(clone@data, inf_df)

    # Extract edge table from PHYLIP output 
    edges <- getPhylipEdges(phylip_out, id_map=id_map)
    
    # Modify PHYLIP tree to remove 0 distance edges
    mod_list <- modifyPhylipEdges(edges, clone)

    # Convert edges and clone data to igraph graph object
    graph <- phylipToGraph(mod_list$edges, mod_list$clone)
    
    return(graph)
}