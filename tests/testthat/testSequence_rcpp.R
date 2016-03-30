file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

test_that("getSeqDistance rcpp reproduces results",{
    # Ungapped examples
    expect_equal(getSeqDistance("ATGGC", "ATGGG", rcpp=T), 1)
    expect_equal(getSeqDistance("ATGGC", "ATG??", rcpp=T), 2)
    
    # Gaps will be treated as Ns with a gap=0 distance matrix
    expect_equal(
        getSeqDistance("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=0), rcpp=T),
        0)
    
    # Gaps will be treated as universally non-matching characters with gap=1
    expect_equal(
        getSeqDistance("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=1), rcpp=T),
        2)
    
    # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
    expect_equal(
        getSeqDistance("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1), rcpp=T),
        1)
    
    # Gaps of equivalent run lengths are not counted as gaps
    expect_equal(
        getSeqDistance("ATG-C", "ATG-C", dist_mat=getDNAMatrix(gap=-1), rcpp=T),
        0)
    
    # Overlapping runs of gap characters are counted as a single gap
    expect_equal(
        getSeqDistance("ATG-C", "AT--C", dist_mat=getDNAMatrix(gap=-1), rcpp=T),
        1) 
    
    expect_equal(
        getSeqDistance("A-GGC", "AT--C", dist_mat=getDNAMatrix(gap=-1), rcpp=T),
        1)
    
    expect_equal(
        getSeqDistance("AT--C", "AT--C", dist_mat=getDNAMatrix(gap=-1), rcpp=T),
        0)
    
    # Discontiguous runs of gap characters each count as separate gaps
    expect_equal(
        getSeqDistance("-TGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1), rcpp=T),
        2)
})

test_that("getSeqMatrix rcpp reproduces restuls", {
    # Gaps will be treated as Ns with a gap=0 distance matrix
    obs <- getSeqMatrix(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
                 dist_mat=getDNAMatrix(gap=0), rcpp=T)
    expect_equal(obs,
                 matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0),ncol=4),
                 check.attributes=F)
    # Gaps will be treated as universally non-matching characters with gap=1
    obs <- getSeqMatrix(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
                 dist_mat=getDNAMatrix(gap=1), rcpp=T)
    expect_equal(obs,
                 matrix(c(0, 1, 1, 2, 1, 0, 0, 3, 1, 0, 0, 3, 2, 3, 3, 0),ncol=4),
                 check.attributes=F)
    
    # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
    obs <- getSeqMatrix(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
                 dist_mat=getDNAMatrix(gap=-1), rcpp=T)
    expect_equal(obs,
                 matrix(c(0, 1, 1, 1, 1, 0, 0, 2, 1, 0, 0, 2, 1, 2, 2, 0),ncol=4),
                 check.attributes=F)
})

test_that("testSeqEqual rcpp reproduces results", {
    # Ignore gaps
    expect_true(testSeqEqual("ATG-C", "AT--C", rcpp=T))
    expect_true(testSeqEqual("ATGGC", "ATGGN", rcpp=T))
    expect_false(testSeqEqual("AT--T", "ATGGC", rcpp=T))
    
    # Ignore only Ns
    expect_false(testSeqEqual("ATG-C", "AT--C", ignore="N", rcpp=T))
    expect_true(testSeqEqual("ATGGC", "ATGGN", ignore="N", rcpp=T))
    expect_false(testSeqEqual("AT--T", "ATGGC", ignore="N", rcpp=T))
})



test_that("rcpp_collapseDuplicates reproduces results", {
    # Example Change-O data.frame
    df <- data.frame(SEQUENCE_ID=LETTERS[1:4],
                     SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
                     TYPE=c("IgM", "IgG", "IgG", "IgA"),
                     SAMPLE=c("S1", "S1", "S2", "S2"),
                     COUNT=1:4,
                     stringsAsFactors=FALSE)
    
    # Annotations are not parsed if neither text_fields nor num_fields is specified
    # The retained sequence annotations will be random
    obs <- rcpp_collapseDuplicates(df, verbose=F)
    exp <- data.frame(
        "SEQUENCE_ID" = c("C", "A"),
        "SEQUENCE_IMGT" = c("NAACTGGN", "CCCCTGGG"),
        "TYPE" = c("IgG","IgM"),
        "SAMPLE" = c("S2", "S1"),
        "COUNT" = c(3,1),
        stringsAsFactors = F
    )
    expect_equal(obs, exp)
    
    # Unique text_fields annotations are combined into a single string with ","
    # num_fields annotations are summed
    # Ambiguous duplicates are discarded
    obs <- rcpp_collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       verbose=F)
    exp$TYPE <- c("IgG","IgG,IgM")
    exp$COUNT <- c(3,3)
    expect_equal(obs, exp)
    
    
    # Use alternate delimiter for collapsing textual annotations
    obs <- rcpp_collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       sep="/", verbose=F)
    exp$TYPE <- c("IgG","IgG/IgM")
    expect_equal(obs, exp)
    
    # Add count of duplicates
    obs <- rcpp_collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       add_count=TRUE, verbose=F)
    exp$TYPE <- c("IgG","IgG,IgM")
    exp$COLLAPSE_COUNT <- c(1,2)
    expect_equal(obs, exp)
    
    # Masking ragged ends may impact duplicate removal
    df$SEQUENCE_IMGT <- maskSeqEnds(df$SEQUENCE_IMGT)
    obs <- rcpp_collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       add_count=TRUE, verbose=F)    
    exp <- data.frame(
        "SEQUENCE_ID" = "A",
        "SEQUENCE_IMGT" = "NNNCTGNN",
        "TYPE" = "IgA,IgG,IgM",
        "SAMPLE" = "S1,S2",
        "COUNT" = 10,
        "COLLAPSE_COUNT" = 4,
        stringsAsFactors = F
    ) 
    expect_equal(obs, exp)
    
})