graph <- ExampleTrees[[5]]

test_that("Test summarizeSubtrees",{
    obs_graph_sum_sub <- summarizeSubtrees(graph, fields="ISOTYPE", root="Germline")
    load("../data-tests/exp_graph_sum_sub.rda")
    expect_equal(obs_graph_sum_sub, exp_graph_sum_sub, tolerance=1)
})

test_that("Test getPathLengths",{
    
    # Consider all nodes
    pl <- getPathLengths(graph, root="Germline")
    expect_equal(pl$STEPS, c(1,2,0,3,2))
    expect_equal(pl$DISTANCE, c(0,1,0,2,1))
    
    # Exclude nodes without an isotype annotation from step count
    graph_2 <- graph
    iso <- V(graph_2)$ISOTYPE
    iso[1] <- NA
    graph_2 <- set_vertex_attr(graph_2, name="ISOTYPE", value=iso)
    pl_exclude <- getPathLengths(graph_2, root="Germline", field="ISOTYPE", 
                                 exclude=NA)
    expect_equal(pl_exclude$STEPS, c(0,1,0,2,1))
    expect_equal(pl_exclude$DISTANCE, c(0,1,0,2,1))
})

test_that("Test getMRCA",{
    mrca <- getMRCA(graph, path="steps", root="Germline")
    expect_equal(mrca$NAME, "GN5SHBT08J0DZP")
    
    mrca_2 <- getMRCA(graph_2, path="distance", root="Germline", 
                      field="ISOTYPE", exclude=NA)
    expect_equal( mrca_2$NAME,  c("GN5SHBT07FM4BU", "GN5SHBT07IBCBZ"))
    
})