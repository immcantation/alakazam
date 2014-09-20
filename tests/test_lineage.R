#' Lineage tests
#' @author  Jason Anthony Vander Heiden  
#' @date    2014.9.20

#### Imports ####
self_path <- dirname(sys.frame(1)$ofile)
test_path <- file.path(dirname(self_path), "R")
source(file.path(test_path, "DataCore.R"), chdir=TRUE)
source(file.path(test_path, "SeqCore.R"), chdir=TRUE)
source(file.path(test_path, "ClonalLineage.R"), chdir=TRUE)

#### Run paramaters ####
clone_file <- "/mnt/data/oconnor_mg_memory/changeo/RQ2341_HD07M/RQ2341_HD07M_db-pass_germ-pass_clone-pass.tab"
text_fields <- "PRCONS"
num_fields <- "DUPCOUNT"

#### Preprocessing ####
data_df <- readChangeoDb(clone_file)
clone_df <- subset(data_df, CLONE == 34)
clone <- prepareClone(clone_df, text_fields=text_fields, num_fields=num_fields)

#### PHYLIP steps ####
temp_path <- makeTempDir(paste0(clone@clone, "-phylip"))
id_map <- writePhylipInput(clone, temp_path)
runPhylip(temp_path)
phylip_out <- readPhylipOutput(temp_path)
checkPhylipOutput(phylip_out)
getPhylipInferred(phylip_out, text_fields=text_fields, num_fields=num_fields)
getPhylipEdges(phylip_out)
# modify edges
# Translate TAXA -> SEQUENCE_ID
# delete temp files
# Return igraph object