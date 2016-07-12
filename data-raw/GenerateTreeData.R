# Generate example trees

# Imports
library(alakazam)
library(dplyr)
library(igraph)

#### Generate graph objects ####

# Load and filter data
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
db <- readChangeoDb(file)

# Preprocess clones
clones <- db %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("SAMPLE", "ISOTYPE"),
                                num_fields="DUPCOUNT"))
# Build lineages
dnapars_exec <- "~/apps/phylip-3.69/dnapars"
graphs <- lapply(clones$CHANGEO, buildPhylipLineage,
                 dnapars_exec=dnapars_exec, rm_temp=TRUE)

# Subset to trees with at least 5 nodes
graphs[sapply(graphs, is.null)] <- NULL
ExampleTrees <- graphs[sapply(graphs, vcount) >= 4]

#### Save to data location ####

devtools::use_data(ExampleTrees, overwrite=TRUE)
