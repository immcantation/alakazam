ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)

test_that("countGenes",{
    # Without copy numbers
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="family")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.76, 0.87, 0.24, 0.05, 0.05, 0.02, 0.01),
                 tolerance=0.001)
    
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="gene")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.41, 0.35, 0.24, 0.33, 0.26, 0.11, 0.10, 0.05, 0.05,
                   0.04, 0.02, 0.02, 0.01, 0.01),
                 tolerance=0.001)
    
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="allele")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.350, 0.325, 0.215, 0.330, 0.085, 0.150, 0.100, 0.100, 0.070,
                   0.050, 0.050, 0.025, 0.040, 0.030, 0.020, 0.020,
                   0.010, 0.010, 0.010, 0.010),
                 tolerance=0.01)
    
    # With copy numbers and multiple groups
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="family")
    expect_equal(genes$SEQ_FREQ, 
                 c(0.78, 0.95, 0.95, 0.22, 0.63, 0.37, 0.60, 0.20, 0.40, 0.67,
                   0.60, 0.80, 0.40, 0.33, 0.03, 0.05, 0.01, 0.01),
                 tolerance=0.01)
    
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="gene")
    expect_equal(round(genes$SEQ_FREQ,2)[1:12], 
                 c(0.61, 0.75, 0.22, 0.54, 0.40, 0.37, 0.17,
                   0.20, 0.15, 0.20, 0.60, 0.20),
                 tolerance=0.01)
    
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="allele")
    expect_equal(genes$SEQ_FREQ[1:12], 
                 c(0.61, 0.72, 0.16, 0.41, 0.40, 0.37, 0.06, 0.20, 0.08, 0.20, 0.50, 0.08),
                 tolerance=0.01)
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