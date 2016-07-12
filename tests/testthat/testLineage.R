ExampleDb <- file.path("..","data-tests","ExampleDb.gz")
db <- readChangeoDb(ExampleDb)

test_that("makeChangeoClone",{
    # Example Change-O data.frame
    db <- data.frame(SEQUENCE_ID=LETTERS[1:4],
                     SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
                     V_CALL="Homsap IGKV1-39*01 F",
                     J_CALL="Homsap IGKJ5*01 F",
                     JUNCTION_LENGTH=2,
                     GERMLINE_IMGT_D_MASK="CCCCAGGG",
                     CLONE=1,
                     TYPE=c("IgM", "IgG", "IgG", "IgA"),
                     COUNT=1:4,
                     stringsAsFactors=FALSE)
    exp <- data.frame(
        "SEQUENCE_ID" = c("C", "A"),
        "SEQUENCE" = c("NAACTGGN", "CCCCTGGG"),
        "TYPE" = c("IgG", "IgG,IgM"),
        "COUNT" = c(3,3),
        "COLLAPSE_COUNT" = c(1,2),
        stringsAsFactors = F
    )
    
    # Without end masking
    clone <- makeChangeoClone(db, text_fields="TYPE", num_fields="COUNT")
    
    expect_true(inherits(clone, "ChangeoClone"))
    expect_equal(clone@clone, "1")
    expect_equal(clone@germline, "CCCCAGGG")
    expect_equal(clone@v_gene, "IGKV1-39")
    expect_equal(clone@j_gene, "IGKJ5")
    expect_equal(clone@junc_len, 2)
    expect_equal(clone@data, exp, tolerance=0.001)
    
    
    # With end masking
    clone <- makeChangeoClone(db, max_mask=3, text_fields="TYPE", num_fields="COUNT")
    exp <- data.frame(
        "SEQUENCE_ID" = "A",
        "SEQUENCE" = c("NNNCTGNN"),
        "TYPE" = c("IgA,IgG,IgM"),
        "COUNT" = c(10),
        "COLLAPSE_COUNT" = c(4),
        stringsAsFactors = F
    )
    
    expect_true(inherits(clone, "ChangeoClone"))
    expect_equal(clone@clone, "1")
    expect_equal(clone@germline, "CCCCAGGG")
    expect_equal(clone@v_gene, "IGKV1-39")
    expect_equal(clone@j_gene, "IGKJ5")
    expect_equal(clone@junc_len, 2)
    expect_equal(clone@data, exp, tolerance=0.001)
    
})

test_that("buildPhylipLineage", {
    # Preprocess clone
    clone <- subset(db, CLONE == 164)
    clone <- makeChangeoClone(clone, text_fields=c("SAMPLE", "ISOTYPE"), 
                              num_fields="DUPCOUNT")
    
    # Run PHYLIP and process output
    
    # Test for error if dnapars_exec doesn't exist
    dnapars_exec <- "~/dummy/phylip-3.69/dnapars"
    if (file.access(dnapars_exec, mode=1) == -1) {
        expect_error(
            graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE),
            "The file ~/dummy/phylip-3.69/dnapars cannot be executed"
        )
    }
    
    dnapars_exec <- Sys.which('dnapars')
    # If dnapars found, run test, else, skip
    if (dnapars_exec!="") {
        graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)
        
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
                     c("name", "sequence", "SAMPLE", "ISOTYPE", "DUPCOUNT", 
                       "COLLAPSE_COUNT", "label"))

        expect_equal(E(graph)$weight,c(1,1,1,0))
        expect_equal(E(graph)$label,c(1,1,1,0))
               
    }
 
})