#ensure older version of sample() used
R_v <- paste(version$major, version$minor,sep=".")
if ( numeric_version(R_v) >= numeric_version("3.6.0") ) {
    RNGkind(sample.kind="Round")   
}

ExampleTrees <- file.path("..", "data-tests", "ExampleTrees.rda")
load(ExampleTrees)

graph <- ExampleTrees[[5]]
graph_2 <- graph

test_that("Test summarizeSubtrees",{
    obs_graph_sum_sub <- summarizeSubtrees(graph, fields="isotype", root="Germline")
    expected <- file.path("..","data-tests", "exp_graph_sum_sub.rda")
    load(expected)
    names(exp_graph_sum_sub)[2] = "isotype"
    expect_equal(obs_graph_sum_sub, exp_graph_sum_sub, tolerance=1)
})

test_that("Test getPathLengths",{
    
    # Consider all nodes
    pl <- getPathLengths(graph, root="Germline")
    expect_equal(pl$STEPS, c(1,0,2,2))
    expect_equal(pl$DISTANCE, c(28,0,32,34))
    
    # Exclude nodes without an isotype annotation from step count
    pl_exclude <- getPathLengths(graph_2, root="Germline", field="isotype", 
                                 exclude=NA)
    expect_equal(pl_exclude$STEPS, c(0,0,1,1))
    expect_equal(pl_exclude$DISTANCE, c(28,0,32,34))
})

test_that("Test getMRCA and testMRCA",{
    mrca <- getMRCA(graph, path="steps", root="Germline")
    expect_equal(mrca$NAME, "Inferred1")
    
    mrca_2 <- getMRCA(graph_2, path="distance", root="Germline", 
                      field="isotype", exclude=NA)
    expect_equal( mrca_2$NAME,  c("GN5SHBT03CABH1"))
    
    graphs <- ExampleTrees[1-10]
     
    # Perform MRCA test on isotypes
    set.seed(8)
    x <- testMRCA(graphs, "isotype", nperm=10)
    x_tests <- slot(x, "tests")
    expect_equal(x_tests$ANNOTATION, c("IgA","IgA,IgG","IgG"))
    expect_equal(x_tests$COUNT, c(16, 1, 31))
    expect_equal(x_tests$EXPECTED, c(13.6, 1.2, 33.8))
    expect_equal(x_tests$PVALUE, c(0, 0.2, 1) )
    
})

test_that("Test tableEdges",{
    graph <- ExampleTrees[[23]]
    
    # Count direct edges between isotypes including inferred nodes
    tab <- tableEdges(graph, "isotype")
    expect_equal(tab$PARENT, c( "IgA", "IgA,IgG", "IgA,IgG", NA))
    expect_equal(tab$CHILD, c( "IgA,IgG", "IgA", "IgG", "IgA"))
    expect_equal(tab$COUNT, c( 1, 1, 3, 1))
    
    # Count direct edges excluding edges to and from germline and inferred nodes
    tab <- tableEdges(graph, "isotype", exclude=c("Germline", NA))
    expect_equal(tab$PARENT, c( "IgA", "IgA,IgG", "IgA,IgG"))
    expect_equal(tab$CHILD, c( "IgA,IgG", "IgA", "IgG"))
    expect_equal(tab$COUNT, c( 1, 1, 3))
    
    # Count indirect edges walking through germline and inferred nodes
    tab <- tableEdges(graph, "isotype", indirect=TRUE, exclude=c("Germline", NA))
    expect_equal(tab$PARENT, c( "IgA", "IgA,IgG", "IgA,IgG"))
    expect_equal(tab$CHILD, c( "IgA,IgG", "IgA", "IgG"))
    expect_equal(tab$COUNT, c( 1, 1, 3))    
         
})

test_that("Test permuteLabels",{
    graph <- ExampleTrees[[23]]
 
    # Permute annotations and plot new tree
    set.seed(43)
    g <- permuteLabels(graph, "isotype")
    expect_equal(V(g)$isotype,
                 c("IgG", "IgG", NA, "IgA", "IgA", "IgG", "IgA,IgG")
    )
})

test_that("Test testEdges", {
    # Define example tree set
    graphs <- ExampleTrees[1-10]
    
    # Perform edge test on isotypes
    set.seed(4)
    x <- slot(testEdges(graphs, "isotype", nperm=10), "tests")
    expect_equal(x$PARENT[1:5], c("IgA", "IgA", "IgA", "IgA,IgG", "IgA,IgG"))
    expect_equal(x$CHILD[5:10], c("IgA,IgG", "IgG", "IgG", "IgA", "IgD,IgG", "IgG"))
    expect_equal(x$COUNT[1:5], c(39, 3, 2, 29, 1))
    expect_equal(x$EXPECTED[5:10], c(2, 6.25, 1.33, 7, 1.00, 133.7), tolerance=0.01)
    expect_equal(x$PVALUE[3:8], c(0.7, 0.0, 0.5, 0.0, 0.0, 1.0), tolerance=0.01)
})
