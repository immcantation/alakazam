# Lineage tests
# @author  Jason Anthony Vander Heiden  
# @date    2014.9.24

#### Run paramaters ####
dnapars_exec <- file.path(Sys.getenv('HOME'), 'apps', 'phylip-3.69', 'dnapars')
clone_file <- "/mnt/data/oconnor_mg_memory/changeo/RQ2341_HD07M/RQ2341_HD07M_db-pass_germ-pass_clone-pass.tab"
text_fields <- "PRCONS"
num_fields <- "DUPCOUNT"

#### Preprocessing ####
data_df <- readChangeoDb(clone_file)
clone_df <- subset(data_df, CLONE == 34)
clone <- prepChangeoClone(clone_df, text_fields=text_fields, num_fields=num_fields)

#### PHYLIP steps ####
graph <- buildPhylipLineage(clone, dnapars_exec=dnapars_exec, rm_temp=TRUE)
ly <- layout.reingold.tilford(graph, root="Germline", circular=F, flip.y=T)
plot(graph, layout=ly)
