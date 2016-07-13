ExampleDb <- system.file("tests/data-tests", "ExampleDb.gz", package="alakazam")
db <- readChangeoDb(ExampleDb)

test_that("calcCoverage", {
    # Calculate clone sizes
    clones <- countClones(db, groups="SAMPLE")
    # Calculate 1st order coverage for a single sample
    obs <- calcCoverage(clones$SEQ_COUNT[clones$SAMPLE == "RL01"])
    expect_equal(obs, 0.1608073, tolerance=0.001)
})

test_that("countClones", {
    # Without copy numbers
    clones <- countClones(db, groups="SAMPLE")
    
    # With copy numbers and multiple groups
    clones <- countClones(db, groups=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")
    
    expect_equal(clones$SEQ_COUNT[1:6], c(1, 2, 1, 1, 2, 1))
    expect_equal(clones$COPY_COUNT[1:6], c(10, 2, 1, 1, 2, 1))
    expect_equal(clones$COPY_FREQ[6:11], rep(0.1, 6), tolerance = 0.001)
})


test_that("estimateAbundance", {
    set.seed(90)
    abund <- estimateAbundance(db, "SAMPLE", nboot=100)
    expect_equal(abund$P[1:6], 
                 c(0.038086, 0.038086, 0.012930, 0.012930, 0.012930, 0.012930),
                 tolerance=0.001)
    expect_equal(abund$LOWER[c(1:3,8:10)],
                 c(0.001102, 0.000786, 0, 0, 0, 0),
                 tolerance = 0.001)
    expect_equal(abund$UPPER[45:50],
                 c(0.00758, 0.00598, 0.00932, 0.00630, 0.00659, 0.00834),
                 tolerance = 0.001)
    expect_equal(abund$RANK[1000:1005], c(36, 37, 38, 39, 40, 41))
})

test_that("calcDiversity", {
    # May define p as clonal member counts
    p <- c(1, 1, 3, 10)
    q <- c(0, 1, 2)
    obs <- calcDiversity(p, q)
    exp <- c(4.000000, 2.594272, 2.027027)
    expect_equal(obs, exp, tolerance = 0.001)
    
    # Or proportional abundance
    p <- c(1/15, 1/15, 1/5, 2/3)
    obs <- calcDiversity(p, q)
    expect_equal(obs, exp, tolerance = 0.001)  
})

test_that("rarefyDiversity", {
    set.seed(5)
    # Group by sample identifier
    div <- rarefyDiversity(db, "SAMPLE", step_q=1, max_q=10, nboot=100)
    obs <- div@data[c(1,3,9,20),]
    exp <- data.frame(
        "GROUP" = c("RL01", "RL01", "RL01", "RL02"),
        "Q" = c(0, 2, 8, 8),
        "D" = c(88.22000, 71.51465, 32.51328, 8.94750),
        "D_SD" = c(3.160936, 8.132310, 10.110378, 2.197165),
        "D_LOWER" = c(82.024680, 55.575615, 12.697307, 4.641136),
        "D_UPPER" = c(94.41532, 87.45369, 52.32926, 13.25386),
        "E" = c(1.0000000, 0.8106399, 0.3685478, 0.1404852),
        "E_LOWER" = c(0.92977420, 0.62996616, 0.14392776, 0.07287072),
        "E_UPPER" = c(1.0702258, 0.9913136, 0.5931678, 0.2080996),
        stringsAsFactors = F
    )
    
    expect_equal(colnames(obs), colnames(exp))
    expect_equal(obs, exp, tolerance=0.001, check.attributes=F)
    
    set.seed <- 25
    # Grouping by isotype rather than sample identifier
    expect_warning(
        div <- rarefyDiversity(db, "ISOTYPE", min_n=40, step_q=1, max_q=10, 
                               nboot=100),
        "Not all groups passed threshold")

    obs <- div@data[c(5,13,19,30),]
    exp <- data.frame(
        "GROUP" = c("IgA", "IgG", "IgG", "IgM"),
        "Q" = c(4, 1, 7, 7),
        "D" = c(10.377177, 7.071540, 2.497792, 40.500567),
        "D_SD" = c(3.1856002, 1.3429885, 0.5335591, 7.1375712),
        "D_LOWER" = c(4.133516, 4.439331, 1.452036, 26.511184),
        "D_UPPER" = c(16.620839, 9.703749, 3.543549, 54.489949),
        "E" = c(0.4008180, 0.5087439, 0.1796973, 0.8343751),
        "E_LOWER" = c(0.1596568, 0.3193763, 0.1044630, 0.5461719),
        "E_UPPER" = c(0.6419791, 0.6981115, 0.2549316, 1.1225783),
        stringsAsFactors = F
    )
    expect_equal(colnames(obs), colnames(exp))
})

test_that("testDiversity", {
    set.seed(3)
    # Groups under the size threshold are excluded and a warning message is issued.
    div <- testDiversity(db, "SAMPLE", q=0, min_n=30, nboot=100)
    expect_equal(div@tests$pvalue, 0)
    expect_equal(div@summary$mean, c(88.10, 63.11), tolerance=0.001)
    
    set.seed(3)
    div <- testDiversity(rbind(db, db), "SAMPLE", q=0, min_n=30, nboot=100)
    expect_equal(div@tests$pvalue, 0.88)
    expect_equal(div@summary$mean, c(78.63, 79.58), tolerance=0.001)
})
