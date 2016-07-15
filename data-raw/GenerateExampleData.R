# Generate example trees

# Imports
library(alakazam)
library(dplyr)
library(igraph)

#### Generate example database ####

# Load data
ExampleDb <- readChangeoDb("data-raw/ExampleDb.gz")
ExampleDb <- ExampleDb[c("SEQUENCE_ID",
                         "SEQUENCE_IMGT",
                         "GERMLINE_IMGT_D_MASK",
                         "V_CALL",
                         "V_CALL_GENOTYPED",
                         "D_CALL",
                         "J_CALL",
                         "JUNCTION",
                         "JUNCTION_LENGTH",
                         "NP1_LENGTH",
                         "NP2_LENGTH",
                         "SAMPLE",
                         "ISOTYPE",
                         "DUPCOUNT",
                         "CLONE")]

# Save
devtools::use_data(ExampleDb, overwrite=TRUE)

#### Generate example trees ####

# Preprocess clones
clones <- ExampleDb %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("SAMPLE", "ISOTYPE"),
                                num_fields="DUPCOUNT", add_count=FALSE))
# Build lineages
dnapars_exec <- "~/apps/phylip-3.69/dnapars"
graphs <- lapply(clones$CHANGEO, buildPhylipLineage,
                 dnapars_exec=dnapars_exec, rm_temp=TRUE)

# Subset to trees with at least 5 nodes
graphs[sapply(graphs, is.null)] <- NULL
ExampleTrees <- graphs[sapply(graphs, vcount) >= 4]

# Save
devtools::use_data(ExampleTrees, overwrite=TRUE)
