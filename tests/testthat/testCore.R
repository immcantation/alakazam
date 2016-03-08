test_that("translateStrings", {
    # Using a vector translation
    strings <- LETTERS[1:5]
    translation <- c("POSITION1"="A", "POSITION5"="E")
    obs <- translateStrings(strings, translation)
    exp <-  c("POSITION1", "B", "C", "D", "POSITION5")
    expect_equal(obs, exp)
    
    # Using a list translation
    translation <- list("1-3"=c("A","B","C"), "4-5"=c("D","E"))
    obs <- translateStrings(strings, translation)
    exp <- c("1-3", "1-3", "1-3", "4-5", "4-5")
    expect_equal(obs, exp)
})


test_that("stoufferMeta", {
    # Define p-value and weight vectors
    p <- c(0.1, 0.05, 0.3)
    w <- c(5, 10, 1)
    
    # Unweighted
    obs <- stoufferMeta(p)
    exp <- c("Z"=1.99232360, "pvalue"=0.02316778)
    expect_equal(obs, exp, tolerance=0.001)
    
    # Weighted
    obs <- stoufferMeta(p, w)    
    exp <- c("Z"=2.08291783, "pvalue"=0.01862936)
    expect_equal(obs, exp, tolerance=0.001)
})
