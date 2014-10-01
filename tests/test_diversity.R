# Diversity tests
# @author  Jason Anthony Vander Heiden  
# @date    2014.9.25

#### Run parameters ####
clone_file <- "inst/extdata/changeo_demo.tab"

#### Preprocessing ####
data_df <- readChangeoDb(clone_file)

#### Diversity steps ####
div_curve <- bootstrapDiversity(data_df, "SAMPLE", nboot=500)
plotDiversityCurve(div_curve)
plotDiversityCurve(div_curve, 
                   axis.text=element_text(size=24), 
                   panel.background=element_rect(fill="black"))

div_test <- testDiversity(data_df, 0, "SAMPLE", nboot=500)
print(div_test)
