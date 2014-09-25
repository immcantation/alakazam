# Diversity tests
# @author  Jason Anthony Vander Heiden  
# @date    2014.9.25

#### Run paramaters ####
clone_file <- "/mnt/data/oconnor_mg_memory/changeo/RQ2341_HD07M/RQ2341_HD07M_db-pass_germ-pass_clone-pass.tab"
text_fields <- "PRCONS"
num_fields <- "DUPCOUNT"

#### Preprocessing ####
data_df <- readChangeoDb(clone_file)

#### PHYLIP steps ####
