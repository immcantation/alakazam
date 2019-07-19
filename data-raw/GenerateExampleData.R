# Generate example trees

# Imports
library(alakazam)
library(dplyr)
library(igraph)

#### Generate example database ####

# Load data
ExampleDb <- readChangeoDb("data-raw/ExampleDb.gz")
ExampleDb <- ExampleDb[c("sequence_id",
                         "sequence_alignment",
                         "germline_alignment_d_mask",
                         "v_call",
                         "v_call_genotyped",
                         "d_call",
                         "j_call",
                         "junction",
                         "junction_length",
                         "np1_length",
                         "np2_length",
                         "sample",
                         "isotype",
                         "duplicate_count",
                         "clone_id")]

# Save
usethis::use_data(ExampleDb, overwrite=TRUE)

#### Generate example trees ####

# Preprocess clones
clones <- ExampleDb %>%
    group_by(clone_id) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("sample", "isotype"),
                                num_fields="duplicate_count", add_count=FALSE))
# Build lineages
dnapars_exec <- "~/apps/phylip-3.69/dnapars"
graphs <- lapply(clones$CHANGEO, buildPhylipLineage,
                 dnapars_exec=dnapars_exec, rm_temp=TRUE)

# Subset to trees with at least 5 nodes
graphs[sapply(graphs, is.null)] <- NULL
airr_trees <- graphs[sapply(graphs, vcount) >= 4]

# Place these trees in the same order as the previously used Change-O trees
load("data-raw/ExampleTreesChangeo.rda")
changeo_trees <- ExampleTrees
changeoClones <- unlist(lapply(ExampleTrees,function(x)x$clone))
airrClones <- unlist(lapply(airr_trees,function(x)x$clone))
ExampleTrees <- airr_trees[order(match(airrClones,changeoClones))]

# Save
usethis::use_data(ExampleTrees, overwrite=TRUE)
