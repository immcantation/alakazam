# Generate example trees

# Imports
library(alakazam)
library(dplyr)
library(igraph)

#### Generate graph objects ####

# Load and filter data
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

# Preprocess clones
clones <- df %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("SAMPLE", "ISOTYPE"),
                                num_fields="DUPCOUNT"))
# Build lineages
dnapars_exec <- "~/apps/phylip-3.69/dnapars"
ExampleTrees <- lapply(clones$CHANGEO, buildPhylipLineage,
                      dnapars_exec=dnapars_exec, rm_temp=TRUE)
ExampleTrees[sapply(ExampleTrees, is.null)] <- NULL

#### Save to data location ####

devtools::use_data(ExampleTrees, overwrite=TRUE)