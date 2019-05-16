ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)

#### calcCoverage ####

test_that("calcCoverage", {
    # Calculate clone sizes
    clones <- countClones(db, group="SAMPLE")
    # Calculate 1st order coverage for a single sample
    obs <- calcCoverage(clones$SEQ_COUNT[clones$SAMPLE == "RL01"])
    expect_equal(obs, 0.1608073, tolerance=0.001)
})

#### countClones ####

test_that("countClones", {
    # Without copy numbers
    clones <- countClones(db, group="SAMPLE")
    expect_equal(clones$SEQ_COUNT[1:6], c(31, 15, 5, 4, 4, 4))
    expect_equal(clones$SEQ_FREQ[1:6], 
                 c(0.15, 0.07, 0.02, 0.04, 0.04, 0.02), tolerance=0.01)
    
    # With copy numbers and multiple groups
    clones <- countClones(db, group=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")
    
    expect_equal(clones$SEQ_COUNT[1:6], c(23, 15, 5, 3, 4, 1))
    expect_equal(clones$COPY_COUNT[1:6], c(53, 43, 24, 11, 11, 10))
    expect_equal(clones$COPY_FREQ[6:11], 
                 c(0.71, 0.05, 0.47, 0.42, 0.04, 0.04),
                 tolerance=0.01)
    
    # Toy database
    db_toy <- tibble::data_frame(SEQUENCE_ID=1:10, 
                                 GROUP=c(rep("A", 3), rep("B", 7)),
                                 CLONE=as.character(c(rep(1, 5), 2, 2, 3, 4, 5)),
                                 COPY=10:1)
    ungrouped_toy <- tibble::data_frame(CLONE=as.character(1:5), 
                                        SEQ_COUNT=as.integer(c(5, 2, 1, 1, 1)),
                                        COPY_COUNT=as.integer(c(sum(10:6), sum(5:4), 3, 2, 1)),
                                        SEQ_FREQ=c(5, 2, 1, 1, 1)/10,
                                        COPY_FREQ=c(sum(10:6), sum(5:4), 3, 2, 1)/sum(10:1))
    grouped_toy <- tibble::data_frame(GROUP=c("A", rep("B", 5)),
                                      CLONE=as.character(c(1, 1:5)), 
                                      SEQ_COUNT=as.integer(c(3, 2, 2, 1, 1, 1)),
                                      COPY_COUNT=as.integer(c(sum(10:8), sum(7:6), sum(5:4), 3, 2, 1)),
                                      SEQ_FREQ=c(3/3, 2/7, 2/7, 1/7, 1/7, 1/7),
                                      COPY_FREQ=c(sum(10:8)/sum(10:8), 
                                                  sum(7:6)/sum(7:1), sum(5:4)/sum(7:1), 3/sum(7:1), 2/sum(7:1), 1/sum(7:1)))
    # Check toy ungrouped
    expect_equal(countClones(db_toy, clone="CLONE", copy="COPY"), ungrouped_toy, tolerance=0.01)
    
    # Check toy grouped
    expect_equal(countClones(db_toy, group="GROUP", clone="CLONE", copy="COPY"), grouped_toy, tolerance=0.01)
})

#### calcInferredDiversity ####

test_that("calcInferredDiversity", {
    # May define p as clonal member counts
    p <- c(1, 1, 3, 10)
    q <- c(0, 1, 2)
    obs <- calcInferredDiversity(p, q)
    exp <- c(5,2.69442045025014,2.05842877124677)
    expect_equal(obs, exp, tolerance = 0.001)
  
})

#### bootstrapDiversity ####

test_that("bootstrapDiversity", {
	
	set.seed(90)
	bootstrap_obj <- estimateAbundance(db, group="SAMPLE", nboot=100)
	expect_equal(bootstrap_obj@abund$P[1:5], 
	             c(0.0026, 0.0025, 0.0016, 0.0014, 9e-04),
	             tolerance=0.001)
	expect_equal(bootstrap_obj@abund$LOWER[c(1:3,8:10)],
	             c(0, 0, 0, 0, 0, 0),
	             tolerance = 0.001)
	expect_equal(bootstrap_obj@abund$UPPER[45:50],
	             c(0.00614, 0.01148, 0.03103, 0.00907, 0.00753, 0.01022),
	             tolerance = 0.001)
	expect_equal(bootstrap_obj@abund$RANK[200:203], c(80, 101, 70, 61))

})

	
#### calculateAlphaDiversity ####

test_that("calculateAlphaDiversity", {	
	set.seed(5)
	bootstrap_obj <- estimateAbundance(db, group="SAMPLE", nboot=100)
	diversity_obj <- calculateAlphaDiversity(bootstrap_obj, step_q=1, max_q=10)
	obs <- diversity_obj@div[c(1,3,9,20),]
    
	exp <- data.frame(
	        "SAMPLE" = c("RL01", "RL01", "RL01", "RL02"),
	        "Q" = c(0, 2, 8, 8),
	        "D_ERROR" = c(107.765026, 1.740357, 1.194441, 1.559210),
	        "D" = c(103.230000, 6.579179, 4.289926, 4.336549),
	        "D_LOWER" = c(0.000000, 4.838821, 3.095485, 2.777339),
	        "D_UPPER" = c(210.995026, 8.319536, 5.484366, 5.895759),
	        "E" = c(1.00000000, 0.06373320, 0.04155697, 0.07837609),
	        "E_LOWER" = c(0.000000, 0.04687418, 0.02998629, 0.05019589),
	        "E_UPPER" = c(2.04393128, 0.08059223, 0.05312764, 0.10655628),
	        stringsAsFactors = F
	    )

	expect_equal(colnames(obs), colnames(exp))

	expect_equal(obs$D, exp$D, tolerance=0.001, check.attributes=F)
	expect_equal(obs$D_ERROR, exp$D_ERROR, tolerance=0.001, check.attributes=F)
	expect_equal(obs$D_LOWER, exp$D_LOWER, tolerance=0.001, check.attributes=F)
	expect_equal(obs$D_UPPER, exp$D_UPPER, tolerance=0.001, check.attributes=F)
	expect_equal(obs$Q, exp$Q, tolerance=0.001, check.attributes=F)
	expect_equal(obs$SAMPLE, exp$SAMPLE, tolerance=0.001, check.attributes=F)
    
    
    set.seed(3)
    bootstrap_obj <- estimateAbundance(db, group="SAMPLE", nboot=100)
    diversity_obj <- calculateAlphaDiversity(bootstrap_obj)
	expect_equal(diversity_obj@test$tests$PVALUE[1:5], 
			c(0.42,0.56,0.76,0.82,0.40), tolerance=0.001)
	expect_equal(diversity_obj@test$summary$MEAN[1:5], 
			c(114.59, 71.86939, 46.32797, 31.92986, 23.940736), tolerance=0.001)
			
})