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


test_that("getSegment", {
    kappa_call <- c("Homsap IGKV1D-39*01 F,Homsap IGKV1-39*02 F,Homsap IGKV1-39*01",
                    "Homsap IGKJ5*01 F")
    
    expect_equal(getAllele(kappa_call),c("IGKV1-39*01", "IGKJ5*01" ))
    
    expect_equal(getAllele(kappa_call, first=FALSE),
                 c("IGKV1-39*01,IGKV1-39*02", "IGKJ5*01"))
    
    expect_equal(getAllele(kappa_call, first=FALSE, strip_d=FALSE),
                 c("IGKV1D-39*01,IGKV1-39*02,IGKV1-39*01", "IGKJ5*01" ))
    
    expect_equal(getGene(kappa_call), c("IGKV1-39", "IGKJ5" ))
    
    expect_equal(getGene(kappa_call, first=FALSE),
                 c("IGKV1-39", "IGKJ5"))
    
    expect_equal(getGene(kappa_call, first=FALSE, strip_d=FALSE),
                 c("IGKV1D-39,IGKV1-39", "IGKJ5"))
    
    expect_equal(getFamily(kappa_call), c("IGKV1", "IGKJ5"))
    
    expect_equal(getFamily(kappa_call, first=FALSE), c("IGKV1", "IGKJ5"))
    
    expect_equal(getFamily(kappa_call, first=FALSE, collapse=FALSE),
                 c("IGKV1,IGKV1,IGKV1", "IGKJ5"))
    
    expect_equal(getFamily(kappa_call, first=FALSE, strip_d=FALSE),
                 c("IGKV1D,IGKV1", "IGKJ5"))
    
    heavy_call <- c("Homsap IGHV1-69*01 F,Homsap IGHV1-69D*01 F", 
                    "Homsap IGHD1-1*01 F", 
                    "Homsap IGHJ1*01 F")
    
    expect_equal(getAllele(heavy_call, first=FALSE),
                 c("IGHV1-69*01", "IGHD1-1*01", "IGHJ1*01" ))
    
    expect_equal(getAllele(heavy_call, first=FALSE, strip_d=FALSE),
                 c("IGHV1-69*01,IGHV1-69D*01", "IGHD1-1*01", "IGHJ1*01"))
    
    expect_equal(getGene(heavy_call, first=FALSE), 
                 c( "IGHV1-69", "IGHD1-1", "IGHJ1"))
    
    expect_equal(getGene(heavy_call, first=FALSE, strip_d=FALSE),
                 c("IGHV1-69,IGHV1-69D", "IGHD1-1", "IGHJ1" ))
})

test_that("sortGenes",{
    genes <- c("IGHV1-69D*01", "IGHV1-69*01", "IGHV4-38-2*01", "IGHV1-69-2*01",
            "IGHV2-5*01", "IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
            "IGHV1-2","IGHV1-2*02","IGHV1-69*02")
    # Sort genes by name
    expected <- c("IGHV1-2", "IGHV1-2*01,IGHV1-2*05", "IGHV1-2*02",
                  "IGHV1-69*01", "IGHV1-69D*01", "IGHV1-69*02", "IGHV1-69-2*01",
                  "IGHV1-NL1*01", "IGHV2-5*01", "IGHV4-38-2*01")
    sorted <- sortGenes(genes)
    expect_equal(sorted, expected)
    
    # Sort genes by position in the locus
    expected_locus <- c("IGHV1-NL1*01", "IGHV1-69-2*01", "IGHV1-69*01", 
                        "IGHV1-69D*01", "IGHV1-69*02", "IGHV4-38-2*01",
                        "IGHV2-5*01", "IGHV1-2", "IGHV1-2*01,IGHV1-2*05", 
                        "IGHV1-2*02")
    sorted_locus <- sortGenes(genes, method="pos")
    expect_equal(sorted_locus, expected_locus)
    
})