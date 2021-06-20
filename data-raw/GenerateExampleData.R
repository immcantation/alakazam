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
                         "germline_alignment",
                         "germline_alignment_d_mask",
                         "rev_comp",
                         "productive",
                         "v_call",
                         "v_call_genotyped",
                         "d_call",
                         "j_call",
                         "c_call",
                         "junction",
                         "junction_length",
                         "np1_length",
                         "np2_length",
                         "duplicate_count",
                         "clone_id",
                         "sample_id")]
c_trans <- c(IGHM="IgM", IGHD="IgD", IGHA="IgA", IGHG="IgG")
ExampleDb <- ExampleDb %>%
    mutate(c_call=translateStrings(c_call, c_trans),
           germline_alignment=germline_alignment_d_mask)

# Save
usethis::use_data(ExampleDb, overwrite=TRUE)

#### Generate example trees ####

# Preprocess clones
clones <- ExampleDb %>%
    group_by(clone_id) %>%
    do(CHANGEO=makeChangeoClone(., id="sequence_id", seq="sequence_alignment", 
                                germ="germline_alignment", text_fields=c("sample_id", "c_call"),
                                num_fields="duplicate_count", add_count=FALSE))
# Build lineages
phylip_exec <- "~/local/apps/phylip-3.695/bin/dnapars"
graphs <- lapply(clones$CHANGEO, buildPhylipLineage,
                 phylip_exec=phylip_exec, rm_temp=TRUE)

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

# Load one sequence db, with all airr fields
SingleDb <- read_rearrangement("data-raw/GN5SHBT02D2WUN_igblast_db-pass.tsv")
usethis::use_data(SingleDb, overwrite=TRUE)
