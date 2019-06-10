ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)

#### calcCoverage ####

test_that("calcCoverage", {
    # Calculate clone sizes
    clones <- countClones(db, groups="SAMPLE")
    # Calculate 1st order coverage for a single sample
    obs <- calcCoverage(clones$SEQ_COUNT[clones$SAMPLE == "RL01"])
    expect_equal(obs, 0.1608073, tolerance=0.001)
})

#### countClones ####

test_that("countClones", {
	# Calculate clone sizes
	clones <- countClones(db, groups="SAMPLE")
	# Calculate 1st order coverage for a single sample
	obs <- calcCoverage(clones$SEQ_COUNT[clones$SAMPLE == "RL01"])
	expect_equal(obs, 0.1608073, tolerance=0.001)

	clones <- countClones(db, groups="SAMPLE")
	expect_equal(clones$SEQ_COUNT[1:6], c(31, 15, 5, 4, 4, 4))
	expect_equal(clones$SEQ_FREQ[1:6], 
	             c(0.15, 0.07, 0.02, 0.04, 0.04, 0.02), tolerance=0.01)

	# With copy numbers and multiple groups
	clones <- countClones(db, groups=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT")

	expect_equal(clones$SEQ_COUNT[1:6], c(23, 15, 5, 3, 4, 1))
	expect_equal(clones$COPY_COUNT[1:6], c(53, 43, 24, 11, 11, 10))
	expect_equal(clones$COPY_FREQ[6:11], 
	             c(0.71, 0.05, 0.47, 0.42, 0.04, 0.04),
	             tolerance=0.01)
    
	# Toy database
	db_toy <- tibble::tibble(SEQUENCE_ID=1:10, 
	                          GROUP=c(rep("A", 3), rep("B", 7)),
	                          CLONE=as.character(c(rep(1, 5), 2, 2, 3, 4, 5)),
	                          COPY=10:1)
	ungrouped_toy <- tibble::tibble(CLONE=as.character(1:5), 
	                                 SEQ_COUNT=as.integer(c(5, 2, 1, 1, 1)),
	                                 COPY_COUNT=as.integer(c(sum(10:6), sum(5:4), 3, 2, 1)),
	                                 SEQ_FREQ=c(5, 2, 1, 1, 1)/10,
	                                 COPY_FREQ=c(sum(10:6), sum(5:4), 3, 2, 1)/sum(10:1))
	grouped_toy <- tibble::tibble(GROUP=c("A", rep("B", 5)),
	                               CLONE=as.character(c(1, 1:5)), 
	                               SEQ_COUNT=as.integer(c(3, 2, 2, 1, 1, 1)),
	                               COPY_COUNT=as.integer(c(sum(10:8), sum(7:6), sum(5:4), 3, 2, 1)),
	                               SEQ_FREQ=c(3/3, 2/7, 2/7, 1/7, 1/7, 1/7),
	                               COPY_FREQ=c(sum(10:8)/sum(10:8), 
	                                           sum(7:6)/sum(7:1), sum(5:4)/sum(7:1), 3/sum(7:1), 2/sum(7:1), 1/sum(7:1)))
	# Check toy ungrouped
	expect_equal(countClones(db_toy, clone="CLONE", copy="COPY"), ungrouped_toy, tolerance=0.01)

	# Check toy grouped
	expect_equal(countClones(db_toy, groups="GROUP", clone="CLONE", copy="COPY"), grouped_toy, tolerance=0.01)
})

#### calcInferredDiversity ####

test_that("calcDiversity", {
	# May define p as clonal member counts
	p <- c(1, 1, 3, 10)
	q <- c(0, 1, 2)
	obs <- calcDiversity(p, q)
	exp <- c(4.000000, 2.594272, 2.027027)
	expect_equal(obs, exp, tolerance=0.001)

	# Or proportional abundance
	p <- c(1/15, 1/15, 1/5, 2/3)
	obs <- calcDiversity(p, q)
	expect_equal(obs, exp, tolerance=0.001)  
})

#### bootstrapDiversity ####

test_that("estimateAbundance", {
	set.seed(90)
	bootstrap_obj <- estimateAbundance(db, group="SAMPLE", nboot=100)
	expect_equal(bootstrap_obj@abundance$P[1:5], 
	             c(0.0372, 0.0370, 0.0139, 0.0133, 0.0126),
	             tolerance=0.001)
	expect_equal(bootstrap_obj@abundance$LOWER[c(1:3,8:10)],
	             c(0.00411, 0, 0, 0, 0, 0),
	             tolerance = 0.001)
	expect_equal(bootstrap_obj@abundance$UPPER[45:50],
	             c(0.00956, 0.00907, 0.00853, 0.00907, 0.00853, 0.00853),
	             tolerance = 0.001)
	expect_equal(bootstrap_obj@abundance$RANK[200:203], c(200, 201, 202, 203))
	
	# Grouping by isotype rather than sample identifier should raise warning
	set.seed(90)
	expect_warning(bootstrap_obj <- estimateAbundance(db, group="ISOTYPE", nboot=100),
	               "Not all groups passed threshold")
})

	
#### alphaDiversity ####

test_that("alphaDiversity", {	
	set.seed(5)
	bootstrap_obj <- estimateAbundance(db, group="SAMPLE", nboot=100)
	diversity_obj <- alphaDiversity(bootstrap_obj, step_q=1, max_q=10)
	obs <- diversity_obj@diversity[c(1,3,9,20),]

	exp <- data.frame(
	        "SAMPLE" = c("RL01", "RL01", "RL01", "RL02"),
	        "Q" = c(0, 2, 8, 8),
	        "D_ERROR" = c(32.96, 13.62, 12.95, 4.73),
	        "D" = c(91.29, 51.66, 26.76, 8.36),
	        "D_LOWER" = c(58.33, 38.04, 13.82, 3.63),
	        "D_UPPER" = c(124.3, 65.3, 39.7, 13.1),
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
	diversity_obj <- alphaDiversity(bootstrap_obj)
	expect_equal(diversity_obj@tests$PVALUE[1:5], 
	        c(0.46,0.48,0.52,0.58,0.74), tolerance=0.001)
	expect_equal(diversity_obj@summary$MEAN[1:5], 
	        c(94, 91.2, 88.5, 85.9, 83.3), tolerance=0.001)
			
})


#### betaDiversity ####

test_that("betaDiversity", {	

	set.seed(3)
	beta_db <- db %>% dplyr::rowwise() %>% dplyr::mutate(RANDOM = sample(1:4, 1, replace=TRUE))
	diversity_obj  <- betaDiversity(beta_db, 
		comparisons=list("1-2"=c("1", "2"), "1-3"=c("1", "3")), group="RANDOM")

	obs <- diversity_obj@diversity[c(1,5,12,72),]

	exp <- data.frame(
	        "COMPARISON" = c("1-2", "1-2", "1-2", "1-3"),
	        "Q" = c(0, 0.4, 1.1, 3.0),
	        "D_ERROR" = c(0.2257, 0.1997, 0.1563, 0.3084),
	        "D" = c(1.790, 1.752, 1.634, 1.275), 
	        "D_LOWER" = c(1.5643, 1.5524, 1.4775, 0.9664),
	        "D_UPPER" = c(2.0157, 1.9517, 1.79017, 1.5831),
	        "E" = c(1.00000000, 0.97879, 0.91273, 0.72458),
	        "E_LOWER" = c(0.87392, 0.867258, 0.82540, 0.54931),
	        "E_UPPER" = c(1.126074, 1.090331, 1.000072, 0.899849),
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
	beta_db <- db %>% dplyr::rowwise() %>% dplyr::mutate(RANDOM = sample(1:4, 1, replace=TRUE))
	diversity_obj  <- betaDiversity(beta_db, comparisons=list("1-2"=c("1", "2"), "1-3"=c("1", "3")), group="RANDOM")
	expect_equal(diversity_obj@tests$PVALUE[1:5], 
	        c(0.5,0.5,0.5,0.5,0.5), tolerance=0.001)
	expect_equal(diversity_obj@summary$MEAN[1:5], 
	        c(1.40, 1.39, 1.38, 1.37, 1.36), tolerance=0.01)
		
})

####
## Reproducibility tests
####

test_that("alphaDiversity reproduces rarefyDiversity", {
    
    # Default params for test in now deprecated rearefyDiversity
    #
    # estimateAbundance(data, group = NULL, clone = "CLONE", copy = NULL,
    #                   ci = 0.95, nboot = 2000, progress = FALSE)
    # rarefyDiversity(data, group, clone = "CLONE", copy = NULL, min_q = 0,
    #                 max_q = 4, step_q = 0.05, min_n = 30, max_n = NULL, ci = 0.95,
    #                 nboot = 2000, uniform = TRUE, progress = FALSE)
    
    # Reproduce old test with alphaDiversity, and expect same results
    set.seed(5)
    diversity_obj <- alphaDiversity(db %>% data.frame, 
                                    group="SAMPLE", clone="CLONE", copy=NULL, 
                                    nboot = 2000,
                                    step_q=0.05, min_q=0, max_q=4, min_n=30, max_n=NULL,
                                    ci=0.95)
    
    obs <- diversity_obj@diversity[c(1,3,9,20),]
    
    # set.seed(5)
    # # Group by sample identifier
    # div <- rarefyDiversity(db, "SAMPLE", step_q=1, max_q=10, nboot=100)
    # obs <- div@data[c(1,3,9,20),]
    
    # expected, from old rarefyDiversity test
    exp <- data.frame(
        "SAMPLE" = c("RL01", "RL01", "RL01", "RL02"),
        "Q" = c(0, 2, 8, 8),
        "D_ERROR" = c(3.160936, 8.132310, 10.110378, 2.197165),
        "D" = c(88.22000, 71.51465, 32.51328, 8.94750),
        "D_LOWER" = c(82.024680, 55.575615, 12.697307, 4.641136),
        "D_UPPER" = c(94.41532, 87.45369, 52.32926, 13.25386),
        "E" = c(1.0000000, 0.8106399, 0.3685478, 0.1404852),
        "E_LOWER" = c(0.92977420, 0.62996616, 0.14392776, 0.07287072),
        "E_UPPER" = c(1.0702258, 0.9913136, 0.5931678, 0.2080996),
        stringsAsFactors = F
    )
    
    expect_equal(colnames(obs), colnames(exp))
    expect_equal(obs %>% data.frame, exp, tolerance=0.001, check.attributes=F)
})
