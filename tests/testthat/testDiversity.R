# example_file <- file.path("tests", "data-tests", "ExampleDb.gz")
example_file <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(example_file)

#### calcCoverage ####

test_that("calcCoverage", {
    # Calculate clone sizes
    clones <- countClones(db, groups="SAMPLE", clone="CLONE")
    # Calculate 1st order coverage for a single sample
    obs <- calcCoverage(clones$seq_count[clones$SAMPLE == "RL01"])
    expect_equal(obs, 0.1608073, tolerance=0.001)
})


#### countClones ####

test_that("countClones", {
	# Calculate clone sizes
	clones <- countClones(db, groups="SAMPLE", clone="CLONE")
	expect_equal(clones$seq_count[1:6], c(31, 15, 5, 4, 4, 4))
	expect_equal(clones$seq_freq[1:6], 
	             c(0.15, 0.07, 0.02, 0.04, 0.04, 0.02), 
	             tolerance=0.01)

	# With copy numbers and multiple groups
	clones <- countClones(db, groups=c("SAMPLE", "ISOTYPE"), copy="DUPCOUNT", clone="CLONE")

	expect_equal(clones$seq_count[1:6], c(23, 15, 5, 3, 4, 1))
	expect_equal(clones$copy_count[1:6], c(53, 43, 24, 11, 11, 10))
	expect_equal(clones$copy_freq[6:11], 
	             c(0.71, 0.05, 0.47, 0.42, 0.04, 0.04),
	             tolerance=0.01)
    
	# Toy database
	db_toy <- tibble::tibble(SEQUENCE_ID=1:10, 
	                          GROUP=c(rep("A", 3), rep("B", 7)),
	                          CLONE=as.character(c(rep(1, 5), 2, 2, 3, 4, 5)),
	                          COPY=10:1)
	ungrouped_toy <- tibble::tibble(CLONE=as.character(1:5), 
	                                 seq_count=as.integer(c(5, 2, 1, 1, 1)),
	                                 copy_count=as.integer(c(sum(10:6), sum(5:4), 3, 2, 1)),
	                                 seq_freq=c(5, 2, 1, 1, 1)/10,
	                                 copy_freq=c(sum(10:6), sum(5:4), 3, 2, 1)/sum(10:1))
	grouped_toy <- tibble::tibble(GROUP=c("A", rep("B", 5)),
	                               CLONE=as.character(c(1, 1:5)), 
	                               seq_count=as.integer(c(3, 2, 2, 1, 1, 1)),
	                               copy_count=as.integer(c(sum(10:8), sum(7:6), sum(5:4), 3, 2, 1)),
	                               seq_freq=c(3/3, 2/7, 2/7, 1/7, 1/7, 1/7),
	                               copy_freq=c(sum(10:8)/sum(10:8), 
	                                           sum(7:6)/sum(7:1), sum(5:4)/sum(7:1), 3/sum(7:1), 2/sum(7:1), 1/sum(7:1)))
	# Check toy ungrouped
	expect_equal(countClones(db_toy, clone="CLONE", copy="COPY"), 
	             ungrouped_toy, 
	             tolerance=0.01)

	# Check toy grouped
	expect_equal(countClones(db_toy, groups="GROUP", clone="CLONE", copy="COPY") %>% ungroup(), 
	             grouped_toy, 
	             tolerance=0.01)
	
	# Test how NAs are handled
	db_some_na <- data.frame(sequence_id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
	                         sample = c("S1", "S2", "S4", "S3", "S2", "S3", "S4", "S4", "S3", "S4", "S0", "S0"),
	                         dupcount = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
	                         clone_id = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4, NA, NA), 
	                         stringsAsFactors=F)
	
	expect_warning(clones <- countClones(db_some_na), NULL)
	expect_equal(clones$seq_freq, 
	             c(0.4, 0.3, 0.2, 0.1),
	             tolerance=0.001)
	
	expect_warning(clones <- countClones(db_some_na, groups = "sample"), NULL)
	expect_equal(clones$seq_count, 
	             c(3, 2, 1, 1, 1, 1, 1),
	             tolerance=0.001)
	
	clones <- countClones(db_some_na, remove_na = FALSE)
	expect_equal(clones$seq_count, 
	             c(4, 3, 2, 2, 1),
	             tolerance=0.001)
	
	clones <- countClones(db_some_na, groups = "sample", remove_na = FALSE)
	expect_equal(clones$seq_count, 
	             c(3, 2, 2, 1, 1, 1, 1, 1),
	             tolerance=0.001)
	
	
	db_all_na <- data.frame(sequence_id = c(1, 2, 3),
	                        sample = c("S1", "S2", "S2"),
	                        dupcount = c(1, 1, 1), 
	                        clone_id = c(NA, NA, NA), 
	                        stringsAsFactors=F)
	
	expect_warning(clones <- countClones(db_all_na), "The column clone_id contains no data")
	expect_equal(nrow(clones), 
	             0,
	             tolerance=0.001)
	
	expect_warning(clones <- countClones(db_all_na, groups = "sample"), "The column clone_id contains no data")
	expect_equal(nrow(clones), 
	             0,
	             tolerance=0.001)
	
	expect_warning(clones <- countClones(db_all_na, remove_na = FALSE), "The column clone_id contains no data")
	expect_equal(clones$seq_count, 
	             3,
	             tolerance=0.001)
	
	expect_warning(clones <- countClones(db_all_na, groups = "sample", remove_na = FALSE), "The column clone_id contains no data")
	expect_equal(clones$seq_count, 
	             c(2, 1),
	             tolerance=0.001)
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


#### estimateAbundance ####

test_that("estimateAbundance-current", {
	set.seed(90)
	abund <- estimateAbundance(db, group="SAMPLE", nboot=100, clone="CLONE")
	expect_equal(abund@abundance$p[1:5], 
	             c(0.0372, 0.0370, 0.0139, 0.0133, 0.0126),
	             tolerance=0.001, check.attributes = FALSE)
	expect_equal(abund@abundance$lower[c(1:3,8:10)],
	             c(0.0041, 0.0004, 0, 0, 0, 0),
	             tolerance = 0.001, check.attributes = FALSE)
	expect_equal(abund@abundance$upper[45:50],
	             c(0.0085, 0.0091, 0.0085, 0.0085, 0.0085, 0.0085),
	             tolerance = 0.001, check.attributes = FALSE)
	expect_equal(abund@abundance$rank[200:203], c(200, 201, 202, 203), check.attributes = FALSE)
	
	# Grouping by isotype rather than sample identifier should raise warning
	set.seed(90)
	expect_warning(abund <- estimateAbundance(db, group="ISOTYPE", nboot=100, clone="CLONE"),
	               "Not all groups passed threshold")
})

	
#### alphaDiversity ####

test_that("alphaDiversity", {
    # Test diversity
    set.seed(5)
	abund <- estimateAbundance(db, group="SAMPLE", nboot=100, clone="CLONE")
	div <- alphaDiversity(abund, step_q=1, max_q=10)
	obs <- data.frame(div@diversity[c(1,3,9,20), ])
	exp <- data.frame("SAMPLE" = c("RL01", "RL01", "RL01", "RL02"),
                      "q" = c(0, 2, 8, 8),
        	          "d" = c(88.4000, 72.0492, 33.8868, 8.9179),
	                  "d_sd" = c(3.4902, 9.3973, 11.9207, 2.1954),
        	          "d_lower" = c(81.5592, 53.6307, 10.5224, 4.6149),
        	          "d_upper" = c(95.2407, 90.4677, 57.2511, 13.2210),
        	          "e" = c(1.0000, 0.8150, 0.3833, 0.1397),
        	          "e_lower" = c(0.9226, 0.6066, 0.1190, 0.0723),
        	          "e_upper" = c(1.0773, 1.0233, 0.6476, 0.2071),
        	          stringsAsFactors = F)
	
	expect_equal(colnames(obs), colnames(exp))
	expect_equal(obs, exp, tolerance=0.001, check.attributes=F)

	# Test diversity p-values
	set.seed(3)
	abund <- estimateAbundance(db, group="SAMPLE", nboot=100, clone="CLONE")
	div <- alphaDiversity(abund, step_q=1, max_q=4)
	expect_equal(div@tests$pvalue[1:5], 
	             c(0, 0, 0, 0, 0), tolerance=0.001)
	expect_equal(div@diversity$d[1:5], 
	             c(88.5600, 82.4607, 72.2387, 60.0029, 50.0501), tolerance=0.001)
	
	# Verify two-steps == one-step
	set.seed(3)
	div_single <- alphaDiversity(db, step_q=1, max_q=4, group="SAMPLE", nboot=100, clone="CLONE")
	expect_equal(div, div_single)
})


#### betaDiversity ####

test_that("betaDiversity", {	
    skip("In development")
	set.seed(3)
	beta_db <- db %>% 
	    dplyr::rowwise() %>% 
	    dplyr::mutate(RANDOM = sample(1:4, 1, replace=TRUE)) %>% 
	    ungroup()
	diversity_obj  <- betaDiversity(beta_db, 
		comparisons=list("1-2"=c("1", "2"), "1-3"=c("1", "3")), group="RANDOM", clone="CLONE")

	obs <- data.frame(diversity_obj@diversity[c(1, 5, 12, 72), ])

	exp <- data.frame(
	        "comparison" = c("1-2", "1-2", "1-2", "1-3"),
	        "q" = c(0, 0.4, 1.1, 3.0),
	        "d" = c(1.790, 1.752, 1.634, 1.275), 
	        "d_sd" = c(0.2257, 0.1997, 0.1563, 0.3084),
	        "d_lower" = c(1.5643, 1.5524, 1.4775, 0.9664),
	        "d_upper" = c(2.0157, 1.9517, 1.79017, 1.5831),
	        "e" = c(1.00000000, 0.97879, 0.91273, 0.72458),
	        "e_lower" = c(0.87392, 0.867258, 0.82540, 0.54931),
	        "e_upper" = c(1.126074, 1.090331, 1.000072, 0.899849),
	        stringsAsFactors = F
	    )

	expect_equal(colnames(obs), colnames(exp))
	expect_equal(obs, exp, tolerance=0.001, check.attributes=F)

	# TODO: fix me. This test isn't right, SAMPLE doen'st exist is exp or obs
	# expect_equal(obs$SAMPLE, exp$SAMPLE, tolerance=0.001, check.attributes=F)

	set.seed(3)
	beta_db <- db %>% 
	    dplyr::rowwise() %>% 
	    dplyr::mutate(RANDOM = sample(1:4, 1, replace=TRUE)) %>% 
	    ungroup()
	diversity_obj  <- betaDiversity(beta_db, comparisons=list("1-2"=c("1", "2"), "1-3"=c("1", "3")), group="RANDOM")
	expect_equal(diversity_obj@tests$pvalue[1:5], 
	        c(0.5,0.5,0.5,0.5,0.5), tolerance=0.001)
	expect_equal(diversity_obj@diversity$d[1:5], 
	        c(1.40, 1.39, 1.38, 1.37, 1.36), tolerance=0.01)
})


#### Reproducibility tests ####

test_that("alphaDiversity reproduces rarefyDiversity and testDiversity", {
    # Default params for test in now deprecated rearefyDiversity
    #
    # estimateAbundance(data, group = NULL, clone = "CLONE", copy = NULL,
    #                   ci = 0.95, nboot = 2000, progress = FALSE)
    # rarefyDiversity(data, group, clone = "CLONE", copy = NULL, min_q = 0,
    #                 max_q = 4, step_q = 0.05, min_n = 30, max_n = NULL, ci = 0.95,
    #                 nboot = 2000, uniform = TRUE, progress = FALSE)
    #
    # The old test
    # set.seed(5)
    # div <- rarefyDiversity(db, "SAMPLE", step_q=1, max_q=10, nboot=100)
    # obs <- div@data[c(1, 3, 9, 20),]
    
    # Reproduce old test with alphaDiversity, and expect same results
    set.seed(5)
    diversity_obj <- alphaDiversity(as.data.frame(db), 
                                    group="SAMPLE", clone="CLONE", copy=NULL, 
                                    step_q=1, min_q=0, max_q=10, min_n=30, max_n=NULL,
                                    ci=0.95, nboot=100)
    obs <- diversity_obj@diversity[c(1, 3, 9, 20), ]
    
    # test the deprecated function
    set.seed(5)
    expect_warning(rarefy_obj <- rarefyDiversity(as.data.frame(db), 
                                                 group="SAMPLE", clone="CLONE", copy=NULL, 
                                                 step_q=1, min_q=0, max_q=10, min_n=30, max_n=NULL,
                                                 ci=0.95, nboot=100),
                  "is deprecated")
    # should match same code run outside the function
    expect_equal(diversity_obj, rarefy_obj)
    
    # Check reprod. old testDiversity
    # Defaults:
    # testDiversity(data, q, group, clone = "CLONE", copy = NULL,
    #               min_n = 30, max_n = NULL, nboot = 2000, progress = FALSE)
    # 
    # The old test:
    # set.seed(3)
    # div <- testDiversity(db, "SAMPLE", q=0, min_n=30, nboot=100)
    # expect_equal(div@tests$PVALUE, 0)
    # expect_equal(div@summary$MEAN, c(88.10, 63.11), tolerance=0.001)

    set.seed(3)
    diversity_obj <- alphaDiversity(as.data.frame(db), 
                                    group="SAMPLE", clone="CLONE", copy=NULL, 
                                    step_q=1, min_q=0, max_q=0, min_n=30, max_n=NULL,
                                    ci=0.95, nboot=100)
    
    # test the deprecated function
    set.seed(3)
    expect_warning(testdiv_obj <- testDiversity(as.data.frame(db), 
                                     group="SAMPLE", clone="CLONE", copy=NULL, 
                                     q=0, min_n=30, max_n=NULL,
                                     ci=0.95, nboot=100),
        "is deprecated")
    
    # should be the same as the code run outside the function
    expect_equal(diversity_obj, testdiv_obj)
})

test_that("estimateAbundance reproduces v0.2.11 results", {
    set.seed(90)
    abund <- estimateAbundance(db, group="SAMPLE", min_n=1, uniform=FALSE, nboot=1000, clone="CLONE")
    # Exact v0.2.11 test results with nboot=100
    # expect_equal(abund@abundance$P[1:6], 
    #              c(0.038086, 0.038086, 0.012930, 0.012930, 0.012930, 0.012930),
    #              tolerance=0.001)
    # expect_equal(abund@abundance$LOWER[c(1:3, 8:10)],
    #              c(0.001102, 0.000786, 0, 0, 0, 0),
    #              tolerance = 0.001)
    # expect_equal(abund@abundance$UPPER[45:50],
    #              c(0.00758, 0.00598, 0.00932, 0.00630, 0.00659, 0.00834),
    #              tolerance = 0.001)
    # v0.2.11 test results with nboot=1000
    expect_equal(abund@abundance$p[1:6],
                 c(0.03808, 0.03808, 0.01293, 0.01293, 0.01293, 0.01293),
                 tolerance=0.005, check.attributes = FALSE)
    expect_equal(abund@abundance$lower[c(1:3, 8:10)],
                 c(0.00049, 0, 0, 0, 0, 0),
                 tolerance = 0.005, check.attributes = FALSE)
    expect_equal(abund@abundance$upper[45:50],
                 c(0.00739, 0.00757, 0.00703, 0.00711,0.00662, 0.00709),
                 tolerance = 0.005, check.attributes = FALSE)
    expect_equal(abund@abundance$rank[1000:1005], c(36, 37, 38, 39, 40, 41))
    
    set.seed(90)
    abund <- estimateAbundance(db[c(1, 289), ], group="SAMPLE", min_n=1, uniform=FALSE, nboot=100, clone="CLONE")
    expect_equal(abund@abundance$lower, c(1, 1), check.attributes = FALSE)
    expect_equal(abund@abundance$upper, c(1, 1), check.attributes = FALSE)
    expect_equal(abund@abundance$rank, c(1, 1), check.attributes = FALSE)
})

test_that("rarefyDiversity reproduces v0.2.11 results", {
    set.seed(5)
    # Group by sample identifier
    expect_warning(div <- rarefyDiversity(db, "SAMPLE", step_q=1, max_q=10, nboot=1000),
                   "is deprecated")
    obs <- data.frame(div@diversity[c(1, 3, 9, 20), ])
    # Exact v0.2.11 test results with nboot=100
    # exp <- data.frame("SAMPLE" = c("RL01", "RL01", "RL01", "RL02"),
    #                   "Q" = c(0, 2, 8, 8),
    #                   "D" = c(88.22000, 71.51465, 32.51328, 8.94750),
    #                   "D_SD" = c(3.160936, 8.132310, 10.110378, 2.197165),
    #                   "D_LOWER" = c(82.024680, 55.575615, 12.697307, 4.641136),
    #                   "D_UPPER" = c(94.41532, 87.45369, 52.32926, 13.25386),
    #                   "E" = c(1.0000000, 0.8106399, 0.3685478, 0.1404852),
    #                   "E_LOWER" = c(0.92977420, 0.62996616, 0.14392776, 0.07287072),
    #                   "E_UPPER" = c(1.0702258, 0.9913136, 0.5931678, 0.2080996),
    #                   stringsAsFactors = FALSE)
    # v0.2.11 test results with nboot=10000
    exp <- data.frame("SAMPLE" = c("RL01", "RL01", "RL01", "RL02"),
                      "q" = c(0, 2, 8, 8),
                      "d" = c(88.16, 71.07, 33.07, 8.83),
                      "d_sd" = c(3.34, 9.78, 11.61, 2.39),
                      "d_lower" = c(81.61, 51.90, 10.32, 4.15),
                      "d_upper" = c(94.72, 90.23, 55.82, 13.53),
                      "e" = c(1.00, 0.81, 0.38, 0.14),
                      "e_lower" = c(0.93, 0.59, 0.12, 0.063),
                      "e_upper" = c(1.07, 1.02, 0.63, 0.21),
                      stringsAsFactors = FALSE)
    
    expect_equal(colnames(obs), colnames(exp))
    expect_equal(obs, exp, tolerance=0.05, check.attributes=F)
    
    # Grouping by isotype rather than sample identifier
    set.seed(25)
    expect_warning(div <- rarefyDiversity(db, "ISOTYPE", min_n=40, step_q=1, max_q=10, nboot=1000),
                   "Not all groups passed threshold")
    obs <- data.frame(div@diversity[c(5, 13, 19, 30), ])
    # Exact v0.2.11 test results with nboot=100
    # exp <- data.frame("ISOTYPE" = c("IgA", "IgG", "IgG", "IgM"),
    #                   "Q" = c(4, 1, 7, 7),
    #                   "D" = c(10.377177, 7.071540, 2.497792, 40.500567),
    #                   "D_SD" = c(3.1856002, 1.3429885, 0.5335591, 7.1375712),
    #                   "D_LOWER" = c(4.133516, 4.439331, 1.452036, 26.511184),
    #                   "D_UPPER" = c(16.620839, 9.703749, 3.543549, 54.489949),
    #                   "E" = c(0.4008180, 0.5087439, 0.1796973, 0.8343751),
    #                   "E_LOWER" = c(0.1596568, 0.3193763, 0.1044630, 0.5461719),
    #                   "E_UPPER" = c(0.6419791, 0.6981115, 0.2549316, 1.1225783),
    #                   stringsAsFactors = FALSE)
    # v0.2.11 test results with nboot=10000
    exp <- data.frame("ISOTYPE" = c("IgA", "IgG", "IgG", "IgM"),
                      "q" = c(4, 1, 7, 7),
                      "d" = c(10.02, 7.26, 2.56, 39.96),
                      "d_sd" = c(2.83, 1.45, 0.49, 6.87),
                      "d_lower" = c(4.47, 4.41, 1.59, 26.49),
                      "d_upper" = c(15.57, 10.10, 3.52, 53.43),
                      "e" = c(0.38, 0.52, 0.18, 0.82),
                      "e_lower" = c(0.17, 0.31, 0.11, 0.55),
                      "e_upper" = c(0.60, 0.72, 0.25, 1.10),
                      stringsAsFactors = FALSE)
    expect_equal(colnames(obs), colnames(exp))
    expect_equal(obs, exp, tolerance=0.05, check.attributes=F)
})

test_that("testDiversity reproduces v0.2.11 results", {
    set.seed(3)
    # Groups under the size threshold are excluded and a warning message is issued.
    expect_warning(div <- testDiversity(db, "SAMPLE", q=0, min_n=30, nboot=1000),
                   "is deprecated")
    # Exact v0.2.11 test results with nboot=100
    # expect_equal(div@tests$PVALUE, 0)
    # expect_equal(div@diversity$D, c(88.10, 63.11), tolerance=0.001)
    # v0.2.11 test results with nboot=10000
    expect_equal(div@tests$pvalue, 0, tolerance=0.05)
    expect_equal(div@diversity$d, c(88.13, 63.57), tolerance=0.05)
    
    set.seed(3)
    expect_warning(div <- testDiversity(rbind(db, db), "SAMPLE", q=0, min_n=30, nboot=1000),
                   "is deprecated")
    # Exact v0.2.11 test results with nboot=100
    # expect_equal(div@tests$PVALUE, 0.88)
    # expect_equal(div@diversity$D, c(78.63, 79.58), tolerance=0.001)
    # v0.2.11 test results with nboot=10000
    expect_equal(div@tests$pvalue, 0.88, tolerance=0.05)
    expect_equal(div@diversity$d, c(78.63, 79.80), tolerance=0.05)
})


#### Duplicated Testing for AIRR FORMAT








# example_file <- file.path("tests", "data-tests", "ExampleDb_airr.gz")
example_file <- file.path("..", "data-tests", "ExampleDb_airr.gz")
db <- readChangeoDb(example_file)

#### calcCoverage ####

test_that("calcCoverage", {
    # Calculate clone sizes
    clones <- countClones(db, groups="sample")
    # Calculate 1st order coverage for a single sample
    obs <- calcCoverage(clones$seq_count[clones$sample == "RL01"])
    expect_equal(obs, 0.1608073, tolerance=0.001)
})


#### countClones ####

test_that("countClones", {
	# Calculate clone sizes
	clones <- countClones(db, groups="sample")
	expect_equal(clones$seq_count[1:6], c(31, 15, 5, 4, 4, 4))
	expect_equal(clones$seq_freq[1:6], 
	             c(0.15, 0.07, 0.02, 0.04, 0.04, 0.02), 
	             tolerance=0.01)

	# With copy numbers and multiple groups
	clones <- countClones(db, groups=c("sample", "isotype"), copy="duplicate_count")

	expect_equal(clones$seq_count[1:6], c(23, 15, 5, 3, 4, 1))
	expect_equal(clones$copy_count[1:6], c(53, 43, 24, 11, 11, 10))
	expect_equal(clones$copy_freq[6:11], 
	             c(0.71, 0.05, 0.47, 0.42, 0.04, 0.04),
	             tolerance=0.01)
    
	# Toy database
	db_toy <- tibble::tibble(SEQUENCE_ID=1:10, 
	                          GROUP=c(rep("A", 3), rep("B", 7)),
	                          CLONE=as.character(c(rep(1, 5), 2, 2, 3, 4, 5)),
	                          COPY=10:1)
	ungrouped_toy <- tibble::tibble(CLONE=as.character(1:5), 
	                                 seq_count=as.integer(c(5, 2, 1, 1, 1)),
	                                 copy_count=as.integer(c(sum(10:6), sum(5:4), 3, 2, 1)),
	                                 seq_freq=c(5, 2, 1, 1, 1)/10,
	                                 copy_freq=c(sum(10:6), sum(5:4), 3, 2, 1)/sum(10:1))
	grouped_toy <- tibble::tibble(GROUP=c("A", rep("B", 5)),
	                               CLONE=as.character(c(1, 1:5)), 
	                               seq_count=as.integer(c(3, 2, 2, 1, 1, 1)),
	                               copy_count=as.integer(c(sum(10:8), sum(7:6), sum(5:4), 3, 2, 1)),
	                               seq_freq=c(3/3, 2/7, 2/7, 1/7, 1/7, 1/7),
	                               copy_freq=c(sum(10:8)/sum(10:8), 
	                                           sum(7:6)/sum(7:1), sum(5:4)/sum(7:1), 3/sum(7:1), 2/sum(7:1), 1/sum(7:1)))
	# Check toy ungrouped
	expect_equal(countClones(db_toy, clone="CLONE", copy="COPY") %>% ungroup(), 
	             ungrouped_toy, 
	             tolerance=0.01)

	# Check toy grouped
	expect_equal(countClones(db_toy, groups="GROUP", clone="CLONE", copy="COPY") %>% ungroup(), 
	             grouped_toy, 
	             tolerance=0.01)
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


#### estimateAbundance ####

test_that("estimateAbundance-current", {
	set.seed(90)
	abund <- estimateAbundance(db, group="sample", nboot=100)
	expect_equal(abund@abundance$p[1:5], 
	             c(0.0372, 0.0370, 0.0139, 0.0133, 0.0126),
	             tolerance=0.001, check.attributes = FALSE)
	expect_equal(abund@abundance$lower[c(1:3,8:10)],
	             c(0.0041, 0.0004, 0, 0, 0, 0),
	             tolerance = 0.001, check.attributes = FALSE)
	expect_equal(abund@abundance$upper[45:50],
	             c(0.0085, 0.0091, 0.0085, 0.0085, 0.0085, 0.0085),
	             tolerance = 0.001, check.attributes = FALSE)
	expect_equal(abund@abundance$rank[200:203], c(200, 201, 202, 203), 
	             check.attributes = FALSE)
	
	# Grouping by isotype rather than sample identifier should raise warning
	set.seed(90)
	expect_warning(abund <- estimateAbundance(db, group="isotype", nboot=100),
	               "Not all groups passed threshold")
})

	
#### alphaDiversity ####

test_that("alphaDiversity", {
    # Test diversity
    set.seed(5)
	abund <- estimateAbundance(db, group="sample", nboot=100)
	div <- alphaDiversity(abund, step_q=1, max_q=10)
	obs <- data.frame(div@diversity[c(1,3,9,20), ])
	exp <- data.frame("sample" = c("RL01", "RL01", "RL01", "RL02"),
                      "q" = c(0, 2, 8, 8),
        	          "d" = c(88.4000, 72.0492, 33.8868, 8.9179),
	                  "d_sd" = c(3.4902, 9.3973, 11.9207, 2.1954),
        	          "d_lower" = c(81.5592, 53.6307, 10.5224, 4.6149),
        	          "d_upper" = c(95.2407, 90.4677, 57.2511, 13.2210),
        	          "e" = c(1.0000, 0.8150, 0.3833, 0.1397),
        	          "e_lower" = c(0.9226, 0.6066, 0.1190, 0.0723),
        	          "e_upper" = c(1.0773, 1.0233, 0.6476, 0.2071),
        	          stringsAsFactors = F)
	
	expect_equal(colnames(obs), colnames(exp))
	expect_equal(obs, exp, tolerance=0.001, check.attributes=F)

	# Test diversity p-values
	set.seed(3)
	abund <- estimateAbundance(db, group="sample", nboot=100)
	div <- alphaDiversity(abund, step_q=1, max_q=4)
	expect_equal(div@tests$pvalue[1:5], 
	             c(0, 0, 0, 0, 0), tolerance=0.001)
	expect_equal(div@diversity$d[1:5], 
	             c(88.5600, 82.4607, 72.2387, 60.0029, 50.0501), tolerance=0.001)
	
	# Verify two-steps == one-step
	set.seed(3)
	div_single <- alphaDiversity(db, step_q=1, max_q=4, group="sample", nboot=100)
	expect_equal(div, div_single)
})


#### betaDiversity ####

test_that("betaDiversity", {	
    skip("In development")
	set.seed(3)
	beta_db <- db %>% 
	    dplyr::rowwise() %>% 
	    dplyr::mutate(RANDOM = sample(1:4, 1, replace=TRUE)) %>% 
	    ungroup()
	diversity_obj  <- betaDiversity(beta_db, 
		comparisons=list("1-2"=c("1", "2"), "1-3"=c("1", "3")), group="RANDOM")

	obs <- data.frame(diversity_obj@diversity[c(1, 5, 12, 72), ])

	exp <- data.frame(
	        "comparison" = c("1-2", "1-2", "1-2", "1-3"),
	        "q" = c(0, 0.4, 1.1, 3.0),
	        "d" = c(1.790, 1.752, 1.634, 1.275), 
	        "d_sd" = c(0.2257, 0.1997, 0.1563, 0.3084),
	        "d_lower" = c(1.5643, 1.5524, 1.4775, 0.9664),
	        "d_upper" = c(2.0157, 1.9517, 1.79017, 1.5831),
	        "e" = c(1.00000000, 0.97879, 0.91273, 0.72458),
	        "e_lower" = c(0.87392, 0.867258, 0.82540, 0.54931),
	        "e_upper" = c(1.126074, 1.090331, 1.000072, 0.899849),
	        stringsAsFactors = F
	    )

	expect_equal(colnames(obs), colnames(exp))
	expect_equal(obs, exp, tolerance=0.001, check.attributes=F)

	# TODO: fix me. This test isn't right, sample doen'st exist is exp or obs
	# expect_equal(obs$sample, exp$sample, tolerance=0.001, check.attributes=F)

	set.seed(3)
	beta_db <- db %>% 
	    dplyr::rowwise() %>% 
	    dplyr::mutate(RANDOM = sample(1:4, 1, replace=TRUE)) %>% 
	    ungroup()
	diversity_obj  <- betaDiversity(beta_db, comparisons=list("1-2"=c("1", "2"), "1-3"=c("1", "3")), group="RANDOM")
	expect_equal(diversity_obj@tests$pvalue[1:5], 
	        c(0.5,0.5,0.5,0.5,0.5), tolerance=0.001)
	expect_equal(diversity_obj@diversity$d[1:5], 
	        c(1.40, 1.39, 1.38, 1.37, 1.36), tolerance=0.01)
})


#### Reproducibility tests ####

test_that("alphaDiversity reproduces rarefyDiversity and testDiversity", {
    # Default params for test in now deprecated rearefyDiversity
    #
    # estimateAbundance(data, group = NULL, clone = "CLONE", copy = NULL,
    #                   ci = 0.95, nboot = 2000, progress = FALSE)
    # rarefyDiversity(data, group, clone = "CLONE", copy = NULL, min_q = 0,
    #                 max_q = 4, step_q = 0.05, min_n = 30, max_n = NULL, ci = 0.95,
    #                 nboot = 2000, uniform = TRUE, progress = FALSE)
    #
    # The old test
    # set.seed(5)
    # div <- rarefyDiversity(db, "sample", step_q=1, max_q=10, nboot=100)
    # obs <- div@data[c(1, 3, 9, 20),]
    
    # Reproduce old test with alphaDiversity, and expect same results
    set.seed(5)
    diversity_obj <- alphaDiversity(as.data.frame(db), 
                                    group="sample", clone="clone_id", copy=NULL, 
                                    step_q=1, min_q=0, max_q=10, min_n=30, max_n=NULL,
                                    ci=0.95, nboot=100)
    obs <- diversity_obj@diversity[c(1, 3, 9, 20), ]
    
    # test the deprecated function
    set.seed(5)
    expect_warning(rarefy_obj <- rarefyDiversity(as.data.frame(db), 
                                                 group="sample", clone="clone_id", copy=NULL, 
                                                 step_q=1, min_q=0, max_q=10, min_n=30, max_n=NULL,
                                                 ci=0.95, nboot=100),
                  "is deprecated")
    # should match same code run outside the function
    expect_equal(diversity_obj, rarefy_obj)
    
    # Check reprod. old testDiversity
    # Defaults:
    # testDiversity(data, q, group, clone = "CLONE", copy = NULL,
    #               min_n = 30, max_n = NULL, nboot = 2000, progress = FALSE)
    # 
    # The old test:
    # set.seed(3)
    # div <- testDiversity(db, "sample", q=0, min_n=30, nboot=100)
    # expect_equal(div@tests$PVALUE, 0)
    # expect_equal(div@summary$MEAN, c(88.10, 63.11), tolerance=0.001)

    set.seed(3)
    diversity_obj <- alphaDiversity(as.data.frame(db), 
                                    group="sample", clone="clone_id", copy=NULL, 
                                    step_q=1, min_q=0, max_q=0, min_n=30, max_n=NULL,
                                    ci=0.95, nboot=100)
    
    # test the deprecated function
    set.seed(3)
    expect_warning(testdiv_obj <- testDiversity(as.data.frame(db), 
                                     group="sample", clone="clone_id", copy=NULL, 
                                     q=0, min_n=30, max_n=NULL,
                                     ci=0.95, nboot=100),
        "is deprecated")
    
    # should be the same as the code run outside the function
    expect_equal(diversity_obj, testdiv_obj)
})

test_that("estimateAbundance reproduces v0.2.11 results", {
    set.seed(90)
    abund <- estimateAbundance(db, group="sample", min_n=1, uniform=FALSE, nboot=1000)
    # Exact v0.2.11 test results with nboot=100
    # expect_equal(abund@abundance$P[1:6], 
    #              c(0.038086, 0.038086, 0.012930, 0.012930, 0.012930, 0.012930),
    #              tolerance=0.001)
    # expect_equal(abund@abundance$LOWER[c(1:3, 8:10)],
    #              c(0.001102, 0.000786, 0, 0, 0, 0),
    #              tolerance = 0.001)
    # expect_equal(abund@abundance$UPPER[45:50],
    #              c(0.00758, 0.00598, 0.00932, 0.00630, 0.00659, 0.00834),
    #              tolerance = 0.001)
    # v0.2.11 test results with nboot=1000
    expect_equal(abund@abundance$p[1:6],
                 c(0.03808, 0.03808, 0.01293, 0.01293, 0.01293, 0.01293),
                 tolerance=0.005, check.attributes = FALSE)
    expect_equal(abund@abundance$lower[c(1:3, 8:10)],
                 c(0.00049, 0, 0, 0, 0, 0),
                 tolerance = 0.005, check.attributes = FALSE)
    expect_equal(abund@abundance$upper[45:50],
                 c(0.00739, 0.00757, 0.00703, 0.00711,0.00662, 0.00709),
                 tolerance = 0.005, check.attributes = FALSE)
    expect_equal(abund@abundance$rank[1000:1005], c(36, 37, 38, 39, 40, 41), 
                 check.attributes = FALSE)
    
    set.seed(90)
    abund <- estimateAbundance(db[c(1, 289), ], group="sample", min_n=1, uniform=FALSE, nboot=100)
    expect_equal(abund@abundance$lower, c(1, 1), check.attributes = FALSE)
    expect_equal(abund@abundance$upper, c(1, 1), check.attributes = FALSE)
    expect_equal(abund@abundance$rank, c(1, 1), check.attributes = FALSE)
})

test_that("rarefyDiversity reproduces v0.2.11 results", {
    set.seed(5)
    # Group by sample identifier
    expect_warning(div <- rarefyDiversity(db, "sample", step_q=1, max_q=10, nboot=1000, clone="clone_id"),
                   "is deprecated")
    obs <- data.frame(div@diversity[c(1, 3, 9, 20), ])
    # Exact v0.2.11 test results with nboot=100
    # exp <- data.frame("sample" = c("RL01", "RL01", "RL01", "RL02"),
    #                   "Q" = c(0, 2, 8, 8),
    #                   "D" = c(88.22000, 71.51465, 32.51328, 8.94750),
    #                   "D_SD" = c(3.160936, 8.132310, 10.110378, 2.197165),
    #                   "D_LOWER" = c(82.024680, 55.575615, 12.697307, 4.641136),
    #                   "D_UPPER" = c(94.41532, 87.45369, 52.32926, 13.25386),
    #                   "E" = c(1.0000000, 0.8106399, 0.3685478, 0.1404852),
    #                   "E_LOWER" = c(0.92977420, 0.62996616, 0.14392776, 0.07287072),
    #                   "E_UPPER" = c(1.0702258, 0.9913136, 0.5931678, 0.2080996),
    #                   stringsAsFactors = FALSE)
    # v0.2.11 test results with nboot=10000
    exp <- data.frame("sample" = c("RL01", "RL01", "RL01", "RL02"),
                      "q" = c(0, 2, 8, 8),
                      "d" = c(88.16, 71.07, 33.07, 8.83),
                      "d_sd" = c(3.34, 9.78, 11.61, 2.39),
                      "d_lower" = c(81.61, 51.90, 10.32, 4.15),
                      "d_upper" = c(94.72, 90.23, 55.82, 13.53),
                      "e" = c(1.00, 0.81, 0.38, 0.14),
                      "e_lower" = c(0.93, 0.59, 0.12, 0.063),
                      "e_upper" = c(1.07, 1.02, 0.63, 0.21),
                      stringsAsFactors = FALSE)
    
    expect_equal(colnames(obs), colnames(exp))
    expect_equal(obs, exp, tolerance=0.05, check.attributes=F)
    
    # Grouping by isotype rather than sample identifier
    set.seed(25)
    expect_warning(div <- rarefyDiversity(db, "isotype", min_n=40, step_q=1, max_q=10, nboot=1000, clone="clone_id"),
                   "Not all groups passed threshold")
    obs <- data.frame(div@diversity[c(5, 13, 19, 30), ])
    # Exact v0.2.11 test results with nboot=100
    # exp <- data.frame("isotype" = c("IgA", "IgG", "IgG", "IgM"),
    #                   "Q" = c(4, 1, 7, 7),
    #                   "D" = c(10.377177, 7.071540, 2.497792, 40.500567),
    #                   "D_SD" = c(3.1856002, 1.3429885, 0.5335591, 7.1375712),
    #                   "D_LOWER" = c(4.133516, 4.439331, 1.452036, 26.511184),
    #                   "D_UPPER" = c(16.620839, 9.703749, 3.543549, 54.489949),
    #                   "E" = c(0.4008180, 0.5087439, 0.1796973, 0.8343751),
    #                   "E_LOWER" = c(0.1596568, 0.3193763, 0.1044630, 0.5461719),
    #                   "E_UPPER" = c(0.6419791, 0.6981115, 0.2549316, 1.1225783),
    #                   stringsAsFactors = FALSE)
    # v0.2.11 test results with nboot=10000
    exp <- data.frame("isotype" = c("IgA", "IgG", "IgG", "IgM"),
                      "q" = c(4, 1, 7, 7),
                      "d" = c(10.02, 7.26, 2.56, 39.96),
                      "d_sd" = c(2.83, 1.45, 0.49, 6.87),
                      "d_lower" = c(4.47, 4.41, 1.59, 26.49),
                      "d_upper" = c(15.57, 10.10, 3.52, 53.43),
                      "e" = c(0.38, 0.52, 0.18, 0.82),
                      "e_lower" = c(0.17, 0.31, 0.11, 0.55),
                      "e_upper" = c(0.60, 0.72, 0.25, 1.10),
                      stringsAsFactors = FALSE)
    expect_equal(colnames(obs), colnames(exp))
    expect_equal(obs, exp, tolerance=0.05, check.attributes=F)
})

test_that("testDiversity reproduces v0.2.11 results", {
    set.seed(3)
    # Groups under the size threshold are excluded and a warning message is issued.
    expect_warning(div <- testDiversity(db, "sample", q=0, min_n=30, nboot=1000, clone="clone_id"),
                   "is deprecated")
    # Exact v0.2.11 test results with nboot=100
    # expect_equal(div@tests$PVALUE, 0)
    # expect_equal(div@diversity$D, c(88.10, 63.11), tolerance=0.001)
    # v0.2.11 test results with nboot=10000
    expect_equal(div@tests$pvalue, 0, tolerance=0.05)
    expect_equal(div@diversity$d, c(88.13, 63.57), tolerance=0.05)
    
    set.seed(3)
    expect_warning(div <- testDiversity(rbind(db, db), "sample", q=0, min_n=30, nboot=1000, clone="clone_id"),
                   "is deprecated")
    # Exact v0.2.11 test results with nboot=100
    # expect_equal(div@tests$PVALUE, 0.88)
    # expect_equal(div@diversity$D, c(78.63, 79.58), tolerance=0.001)
    # v0.2.11 test results with nboot=10000
    expect_equal(div@tests$pvalue, 0.88, tolerance=0.05)
    expect_equal(div@diversity$d, c(78.63, 79.80), tolerance=0.05)
})

