file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)
test_that("countGenes",{
    # Without copy numbers
    genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="family")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.87, 0.05, 0.05, 0.02, 0.01, 0.76, 0.24),
                 tolerance=0.001)
    
    genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="gene")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.33, 0.26, 0.11, 0.10, 0.05, 0.05, 0.04, 0.02, 0.02, 
                   0.01, 0.01, 0.41, 0.35, 0.24),
                 tolerance=0.001)
    
    genes <- countGenes(df, gene="V_CALL", groups="SAMPLE", mode="allele")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.33, 0.15, 0.1, 0.1, 0.07, 0.05, 0.05, 0.04, 0.03, 0.02, 
                   0.02, 0.01, 0.01, 0.01, 0.01, 0.35, 0.325, 0.215, 0.085, 0.025),
                 tolerance=0.001)
    
    # With copy numbers and multiple groups
    genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="family")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.2, 0.8, 0.6, 0.4, 0.6, 0.4, 0.947, 0.027, 0.013, 0.013,
                   0.783, 0.217, 0.667, 0.333, 0.95, 0.05, 0.632, 0.368),
                 tolerance=0.001)
    
    genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="gene")
    expect_equal(genes$SEQ_FREQ[1:12], 
                 c(0.2, 0.6, 0.2, 0.5, 0.4, 0.1, 0.6, 0.4, 0.4, 0.2, 0.1466, 0.12),
                 tolerance=0.001)
    
    genes <- countGenes(df, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="allele")
    expect_equal(genes$SEQ_FREQ[1:12], 
                 c(0.2, 0.6, 0.2, 0.4, 0.3, 0.2, 0.1, 0.5, 0.4, 0.1, 0.4, 0.12),
                 tolerance=0.001)
})