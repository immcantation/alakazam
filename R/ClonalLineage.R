#' Ig lineage reconstruction via maximum parsimony
#' 
#' @author     Jason Anthony Vander Heiden
#' @copyright  Copyright 2014 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @version    0.2.0
#' @date       2014.9.16

#### TODO ####
# Multifurcating tree
#ape::multi2di
# Read tree file into phylo object
#ape::read.tree
# Make/score trees with parsimony
#phangorn::parsimony

#### Imports ####
#self_path <- dirname(sys.frame(1)$ofile)
#source(file.path(self_path, "DataCore.R"), chdir=TRUE)
#source(file.path(self_path, "SeqCore.R"), chdir=TRUE)


#### Constants ####
DNAPARS_EXEC <- file.path(Sys.getenv('HOME'), 'apps', 'phylip-3.69', 'dnapars')


#### Class definitions ####

#' S4 class defining a clone
setClass("ChangeoClone", contains="data.frame",
         slots=c(clone="character", 
                 germline="character", 
                 vgene="character", 
                 jgene="character", 
                 junc_len="numeric"))


#### Preprocessing functions ####

#' Performs preprocessing of a clonal group for tree construction
#' by masking gap positions, masking ragged ends, and removing duplicates
#'
#' @param   data         a data.frame containing the Change-O data for a clone
#' @param   max_mask     the maximum number of characters to mask from the ends
#'                       if NULL no threshold is set
#' @param   text_fields  additional text annotation columns to process during collapse
#' @param   num_fields   additional numeric annotation columns to process during collapse
#' @return  a ChangeoClone object containing the modified clone
prepareClone <- function(data, max_mask=NULL, text_fields=NULL, num_fields=NULL) {
    # Replace gaps with Ns and masked ragged ends
    tmp_df <- data[, c("SEQUENCE_ID", "SEQUENCE_GAP", text_fields, num_fields)]
    tmp_df[, "SEQUENCE_GAP"] <- maskSeqGaps(tmp_df[, "SEQUENCE_GAP"], outer_only=FALSE)
    tmp_df[, "SEQUENCE_GAP"] <- maskSeqEnds(tmp_df[, "SEQUENCE_GAP"], max_mask=max_mask, 
                                            trim=FALSE)
    
    # Remove duplicates
    tmp_df <- collapseDuplicates(tmp_df, id="SEQUENCE_ID", seq="SEQUENCE_GAP",
                                 text_fields=text_fields, num_fields=num_fields)
    
    # Define return object
    clone <- new("ChangeoClone", tmp_df, 
                 clone=as.character(data[1, "CLONE"]),
                 germline=maskSeqGaps(data[1, "GERMLINE_GAP_D_MASK"], outer_only=FALSE), 
                 vgene=getGene(data[1, "V_CALL"]), 
                 jgene=getGene(data[1, "J_CALL"]), 
                 junc_len=data[1, "JUNCTION_GAP_LENGTH"])
    
    return(clone)
}


#### PHYLIP functions ####

#' Create PHYLIP input files in a temporary folder
#'
#' @param   clone  a ChangeoClone object
#' @param   path   a directory to store the write the output files to
#' @return  a named vector translating sequence IDs (names) to PHYLIP taxa (values)
writePhylipInput <- function(clone, path) {
    # Define PHYLIP columns
    v1 <- c(sprintf('%-9s', nrow(clone)),
            sprintf("%-9s", "GERMLINE"), 
            sprintf("SAM%-6s", 1:nrow(clone)))
    v2 <- c(nchar(clone@germline),
            clone@germline, 
            clone[, "SEQUENCE_GAP"])
    phy_df <- data.frame(v1, v2, stringsAsFactors=F)
    
    # Define names vector mapping taxa names to original sequence identifiers
    id_map <- setNames(str_trim(v1[-(1:2)]), clone[, "SEQUENCE_ID"])
    
    # Create PHYLIP input file
    write.table(phy_df, file=file.path(path, "infile"), 
                quote=F, sep=" ", col.names=F, row.names=F)    
    
    return(id_map)
}


#' Run PHYLIP dnapars application
#'
#' @param   path          the temporary directory containing infile
#' @param   dnapars_exec  the path to the dnapars executable
#' @return  NULL
runPhylip <- function(path, dnapars_exec=DNAPARS_EXEC) {
    # Remove old files
    if (file.exists(file.path(path, "outfile"))) { file.remove(file.path(path, "outfile")) }
    if (file.exists(file.path(path, "outtree"))) { file.remove(file.path(path, "outtree")) }    

    # Run dnapars
    phy_options <- c("S", "Y", "I", "4", "5", ".")
    system(paste0("cd ", path, "; ", dnapars_exec), input=c(phy_options, "Y"))
}


#' Reads in the PHYLIP outfile
#'
#' @param   path  the temporary folder containing the dnapars outfile
#' @return  a character vector with each item as a line in the outfile
readPhylipOutput <- function(path) {
    phylip_out <- scan(file.path(path, "outfile"), what="character", sep="\n", 
                       blank.lines.skip=F, strip.white=F)
    return(phylip_out)
}


#' Test for successful PHYLIP dnapars run by checking the outfile
#'
#' @param   phylip_out  a character vector returned by readPhylipOut
#' @return  TRUE if trees built 
#'          FALSE if no trees built
checkPhylipOutput <- function(phylip_out) {
    # Check for failed tree build
    result <- !(any(grepl('-1 trees in all found', phylip_out)))
    
    return(result)
}


#' Extracts inferred sequences from PHYLIP dnapars outfile
#'
#' @param   phylip_out   a character vector returned by readPhylipOutput
#' @param   text_fields  a vector of text columns to fill with default values 
#' @param   num_fields   a vector of numeric columns to fill with default values 
#' @return  data.frame of inferred sequences with columns (TAXA, SEQUENCE)
getPhylipInferred <- function(phylip_out, text_fields=NULL, num_fields=NULL) {
    # Process dnapars output
    seq_start <- min(grep("From\\s+To\\s+Any Steps\\?\\s+State at upper node", 
                          phylip_out, perl=T, fixed=F))
    seq_empty <- grep("^\\s*$", phylip_out[seq_start:length(phylip_out)], perl=T, fixed=F)
    seq_len <- seq_empty[min(which(seq_empty[-1] == (seq_empty[-length(seq_empty)] + 1)))]
    seq_block <- paste(phylip_out[(seq_start + 2):(seq_start + seq_len - 2)], collapse="\n")
    seq_df <- read.table(textConnection(seq_block), as.is=T, fill=T, blank.lines.skip=F)
    
    # Correct first line of block and remove blank rows
    fix.row <- c(1, which(is.na(seq_df[,1])) + 1)
    seq_df[fix.row, ] <- cbind(0, seq_df[fix.row, 1], "no", seq_df[fix.row, 2:5], stringsAsFactors=F)
    seq_df <- seq_df[-(fix.row[-1] - 1), ]
    
    # Create data.frame of inferred sequences
    inferred_names <- unique(grep("^[0-9]+$", seq_df[, 2], value=T))
    inferred_seq <- sapply(inferred_names, function(n) { paste(t(as.matrix(seq_df[seq_df[, 2] == n, -c(1:3)])), collapse="") })
    inferred_df <- data.frame(TAXA=inferred_names, 
                              SEQUENCE=inferred_seq, 
                              stringsAsFactors=F)
    
    # Add additional fields
    if (!is.null(text_fields)) { inferred_df[, text_fields] <- "Inferred" }
    if (!is.null(num_fields)) { inferred_df[, num_fields] <- 0 }
    
    return(inferred_df)
}


#' Extracts graph edge list from a PHYLIP dnapars outfile
#'
#' @param    phylip_out  a character vector returned by readPhylipOutput
#' @returns  a data.frame of edges with columns (from, to, weight)
getPhylipEdges <- function(phylip_out) {
    # Process dnapars output
    edge_start <- min(grep('between\\s+and\\s+length', phylip_out, 
                           perl=TRUE, fixed=FALSE))
    edge_len <- min(grep('^\\s*$', phylip_out[edge_start:length(phylip_out)], 
                         perl=TRUE, fixed=FALSE))
    edge_block <- paste(phylip_out[(edge_start + 2):(edge_start + edge_len - 2)], collapse='\n')
    edge_df <- read.table(textConnection(edge_block), col.names=c('from', 'to', 'weight'), 
                          as.is=TRUE)
    
    return(edge_df)
}


#' Modify edges of phylip output
#'
#' @param   edges    data.frame of edges returned by getPhylipEdges
#' @param   clone    a ChangeoClone object containg sequence data
#' @param   nuc_mat  nucleotide character distance matrix
#' @return  modified edges data.frame
modifyPhylipEdges <- function(edges, clone, nuc_mat=getNucMatrix()) {
    # Move germline to root position
    germline.idx <- which(edges$to == "GERMLINE")
    edges[germline.idx, c('from', 'to')] <- edges[germline.idx, c('to', 'from')]
    
    # Calculate edge mutations
    for (i in 1:nrow(edges)) {
        if (edges$from[i] == "GERMLINE") {
            seq1 <- clone@germline
        } else {
            seq1 <- clone[clone[, "SEQUENCE_ID"] == edges$from[i], "SEQUENCE_GAP"]
        }
        seq2 <- clone[clone[, "SEQUENCE_ID"] == edges$to[i], "SEQUENCE_GAP"]
        edges$weight[i] <- getSeqDistance(seq1, seq2, nuc_mat)        
    }
    
    # Find rows zero weight edges with inferred parent nodes
    remove.row <- which(edges$weight == 0 & 
                        edges$from != "GERMLINE" & 
                        grepl('^\\d+$', edges$from))
    
    # Replace inferred parent nodes with child nodes when edge weight is zero
    while (length(remove.row) > 0) {
        # Remove first node with zero distance to parent
        r <- remove.row[1]
        r.idx <- which(edges[c('from', 'to')] == edges$from[r], arr.ind=T)
        edges[r.idx] <- edges$to[r]
        
        # Recalculate edge weights for modified rows
        r.mod <- r.idx[, 1][r.idx[, 1] != r]
        for (i in r.mod) {
            if (edges$from[i] == "GERMLINE") {
                seq1 <- clone@germline
            } else {
                seq1 <- clone[clone[, "SEQUENCE_ID"] == edges$from[i], "SEQUENCE_GAP"]
            }
            seq2 <- clone[clone[, "SEQUENCE_ID"] == edges$to[i], "SEQUENCE_GAP"]
            edges$weight[i] <- getSeqDistance(seq1, seq2, nuc_mat)      
        }
        
        # Remove row
        edges <- edges[-r, ]
        
        # Re-determine rows to remove
        remove.row <- which(edges$weight == 0 & 
                            edges$from != "GERMLINE" & 
                            grepl('^\\d+$', edges$from))         
    }
    
    return(edges)
}

#' Convert edge data.frame and clone object to igraph graph object
#'
#' @param   edges    data.frame of edges returned by getPhylipEdges
#' @param   clone    a ChangeoClone object containg sequence data
#' @return  an igraph graph object
convertPhylip <- function(edges, clone) {
    # Create igraph object
    g <- graph.data.frame(edges, directed=T)
    V(g)$number <- match(V(g)$name, clone[, "SEQUENCE_ID"])
    c_df <- clone[V(g)$number, ]
}