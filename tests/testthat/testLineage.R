# ExampleDb <- file.path("tests", "data-tests", "ExampleDb_airr.gz")
ExampleDb <- file.path("..", "data-tests", "ExampleDb_airr.gz")
db <- readChangeoDb(ExampleDb)

test_that("makeChangeoClone",{
    # Example AIRR data.frame
    db <- data.frame(sequence_id=LETTERS[1:4],
                     sequence_alignment=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
                     v_call="Homsap IGKV1-39*01 F",
                     j_call="Homsap IGKJ5*01 F",
                     junction_length=2,
                     germline_alignment="CCCCAGGG",
                     clone_id=1,
                     TYPE=c("IgM", "IgG", "IgG", "IgA"),
                     COUNT=1:4,
                     stringsAsFactors=FALSE)
    
    # Example Change-O data.frame
   changeodb <- data.frame(SEQUENCE_ID=LETTERS[1:4],
                  SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
                  V_CALL="Homsap IGKV1-39*01 F",
                  J_CALL="Homsap IGKJ5*01 F",
                  JUNCTION_LENGTH=2,
                  GERMLINE_IMGT_D_MASK="CCCCAGGG",
                  CLONE=1,
                  TYPE=c("IgM", "IgG", "IgG", "IgA"),
                  COUNT=1:4,
                  stringsAsFactors=FALSE)
                   
    exp <- data.frame("sequence_id"=c("C", "A"),
                      "sequence"=c("NAACTGGN", "CCCCTGGG"),
                      "TYPE"=c("IgG", "IgG,IgM"),
                      "COUNT"=c(3, 3),
                      "collapse_count"=c(1, 2),
                      stringsAsFactors=FALSE)
    
    # Without end masking
    clone <- makeChangeoClone(db, text_fields="TYPE", num_fields="COUNT")
    
    changeoclone <- makeChangeoClone(changeodb, id="SEQUENCE_ID", seq="SEQUENCE_IMGT", 
           germ="GERMLINE_IMGT_D_MASK", v_call="V_CALL", j_call="J_CALL", 
           junc_len="JUNCTION_LENGTH", clone="CLONE", text_fields="TYPE", num_fields="COUNT")

    expect_true(inherits(clone, "ChangeoClone"))
    expect_equal(clone@clone, "1")
    expect_equal(clone@germline, "CCCCAGGG")
    expect_equal(clone@v_gene, "IGKV1-39")
    expect_equal(clone@j_gene, "IGKJ5")
    expect_equal(clone@junc_len, 2)
    expect_equal(clone@data, exp, tolerance=0.001)
    expect_equal(changeoclone@data, exp, tolerance=0.001)
    
    # With padding
    db_trim <- db
    db_trim$sequence_alignment <- c("CCCCTGGG", "CCCCTGG", "NAACTGG", "NNNCTG")
    clone <- makeChangeoClone(db_trim, text_fields="TYPE", num_fields="COUNT", pad_end=TRUE)

    db_trim <- changeodb
    db_trim$sequence_alignment <- c("CCCCTGGG", "CCCCTGG", "NAACTGG", "NNNCTG")
    changeoclone <- makeChangeoClone(db_trim, id="SEQUENCE_ID", seq="SEQUENCE_IMGT", 
           germ="GERMLINE_IMGT_D_MASK", v_call="V_CALL", j_call="J_CALL", 
           junc_len="JUNCTION_LENGTH", clone="CLONE",text_fields="TYPE", 
           num_fields="COUNT", pad_end=TRUE)

    expect_true(inherits(clone, "ChangeoClone"))
    expect_equal(clone@clone, "1")
    expect_equal(clone@germline, "CCCCAGGG")
    expect_equal(clone@v_gene, "IGKV1-39")
    expect_equal(clone@j_gene, "IGKJ5")
    expect_equal(clone@junc_len, 2)
    expect_equal(clone@data, exp, tolerance=0.001)
    expect_equal(changeoclone@data, exp, tolerance=0.001)

    # With end masking
    clone <- makeChangeoClone(db, max_mask=3, text_fields="TYPE", num_fields="COUNT")
    changeoclone <- makeChangeoClone(changeodb, max_mask=3, id="SEQUENCE_ID",
            seq="SEQUENCE_IMGT", germ="GERMLINE_IMGT_D_MASK", v_call="V_CALL", 
            j_call="J_CALL", junc_len="JUNCTION_LENGTH", clone="CLONE", 
            text_fields="TYPE", num_fields="COUNT")

    exp <- data.frame("sequence_id"="A",
                      "sequence"=c("NNNCTGNN"),
                      "TYPE"=c("IgA,IgG,IgM"),
                      "COUNT"=c(10),
                      "collapse_count"=c(4),
                      stringsAsFactors=F)
    
    expect_true(inherits(clone, "ChangeoClone"))
    expect_equal(clone@clone, "1")
    expect_equal(clone@germline, "CCCCAGGG")
    expect_equal(clone@v_gene, "IGKV1-39")
    expect_equal(clone@j_gene, "IGKJ5")
    expect_equal(clone@junc_len, 2)
    expect_equal(clone@data, exp, tolerance=0.001)
    expect_equal(changeoclone@data, exp, tolerance=0.001)
})

test_that("buildPhylipLineage", {
    # Preprocess clone
    clone <- subset(db, clone_id == 164)
    clone <- makeChangeoClone(clone, germ="germline_alignment_d_mask", 
                              text_fields=c("sample", "isotype"), 
                              num_fields="duplicate_count")
    
    # Run PHYLIP and process output
    
    # Test for error if phylip_exec doesn't exist
    phylip_exec <- "~/dummy/phylip-3.69/dnapars"
    if (file.access(phylip_exec, mode=1) == -1) {
        expect_error(
            graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE),
            "The file ~/dummy/phylip-3.69/dnapars cannot be executed"
        )
    }
    
    phylip_exec <- Sys.which('dnapars')
    # If dnapars found, run test, else, skip
    if (phylip_exec!="") {
        graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
        
        expect_true(inherits(graph, "igraph"))
        expect_equal(igraph::vcount(graph), 5)
        expect_equal(igraph::ecount(graph), 4)
        expect_true(igraph::is.directed(graph))
        
        expect_equal(igraph::graph_attr_names(graph),
                     c("clone", "v_gene", "j_gene", "junc_len"))
        
        expect_equal(igraph::graph_attr(graph),
                     list("clone"="164", "v_gene"="IGHV3-48","j_gene"="IGHJ2",
                       "junc_len"=66))
        
        expect_equal(V(graph)$name,
                     c("GN5SHBT08J0DZP", "GN5SHBT07FM4BU", "Germline",
                       "GN5SHBT07IAWQ9", "GN5SHBT07IBCBZ"))
        expect_equal(V(graph)$name, V(graph)$label)
        expect_equal(V(graph)$sequence,
                     c("GAGGTGCAGCTGGTGGAGTCTGGGGGANNNGGCTTGGTACAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCNNNNNNNNNNNNAGTAGTTATGAAATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGTNNNNNNGGTAGTACCATATACTACGCAGACTCTGTGAAGNNNGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGTGCGAGAGATTTGGCGCTCGATAGTAGTGGTTATTACCCTAGCTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG",
                       "GAGGTGCAGCTGGTGGAGTCTGGGGGANNNGGCTTGGTACAGCCTGGAGGGTCCCCGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCNNNNNNNNNNNNAGTAGTTATGAAATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGTNNNNNNGGTAGTACCATATACTACGCAGACTCTGTGAAGNNNGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGTGCGAGAGATTTGGCGCTCGATAGTAGTGGTTATTACCCTAGCTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG",
                       "GAGGTGCAGCTGGTGGAGTCTGGGGGANNNGGCTTGGTACAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCNNNNNNNNNNNNAGTAGTTATGAAATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGTNNNNNNGGTAGTACCATATACTACGCAGACTCTGTGAAGNNNGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGTGCGAGAGANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG",
                       "GAGGTGCAGCTGGTGGAGTCTGGGGGANNNGGCTTGGTACAGCCTGGAGGGTCCCCGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCNNNNNNNNNNNNAGTAGTTATGAAATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGGTGGAGTGGGTTTCATACATTAGTAGTAGTNNNNNNGGTAGTACCATATACTACGCAGACTCTGTGAAGNNNGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGTGCGAGAGATTTGGCGCTCGATAGTAGTGGTTATTACCCTAGCTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG",
                       "GAGGTGCAGCTGGTGGAGTCTGGGGGANNNGGCTTGGTACAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCNNNNNNNNNNNNAGTAGTTATGAAATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGTNNNNNNGGTAGTACCATATACTACGCAGACTCTGTGAAGNNNGGCCGATCCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGTGCGAGAGATTTGGCGCTCGATAGTAGTGGTTATTACCCTAGCTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG"
                     ))
        
        expect_equal(igraph::vertex_attr_names(graph),
                     c("name", "sequence", "sample", "isotype", "duplicate_count", 
                       "collapse_count", "label"))

        expect_equal(E(graph)$weight,c(1,1,1,0))
        expect_equal(E(graph)$label,c(1,1,1,0))
               
    }
 
})
