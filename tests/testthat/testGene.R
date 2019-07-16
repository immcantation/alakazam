ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)

### countGenes ####

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

    #Testing of fill

    # Without copy numbers
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="family", fill = T)
    expect_equal(genes$SEQ_FREQ, 
                 c(0.02,0.87,0.05,0.05,0.01,0,0.76,0,0.24,0),
                 tolerance=0.001)

    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="gene", fill = T)
    expect_equal(genes$SEQ_FREQ, 
                 c(0.02,0.04,0.01,0.1,0.11,0.26,0.33,0.02,0,0.05,0.05,0.01,0,0,
                    0,0,0.41,0,0,0,0.35,0,0.24,0),
                 tolerance=0.001)

    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="allele", fill = T)
    expect_equal(genes$SEQ_FREQ, 
                 c(0.02,0.01,0.03,0.01,0.1,0.07,0.04,0.01,0.1,0.15,0.33,0.02,0,
                    0.05,0.05,0,0.01,0,0,0,0,0,0.085,0.325,0,0,0,0,0,0.35,0,0.215,0.025,0),
                 tolerance=0.01)

    # With copy numbers and multiple groups
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="family", fill = T)
    expect_equal(round(genes$SEQ_FREQ,2)[1:12], 
                 c(0,0.8,0,0.2,0,0,0.6,0.4,0,0,0,0.6),
                 tolerance=0.01)

    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="gene", fill = T)
    expect_equal(round(genes$SEQ_FREQ,2)[1:12], 
                 c(0,0,0,0.2,0,0,0.6,0,0,0,0.2,0),
                 tolerance=0.01)

    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="allele", fill = T)
    expect_equal(genes$SEQ_FREQ[1:12], 
                 c(0,0,0,0,0.2,0,0,0,0,0,0.6,0),
                 tolerance=0.01)
})

### getSegment ####

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
                 c("IGHV1-69", "IGHD1-1", "IGHJ1"))
    
    expect_equal(getGene(heavy_call, first=FALSE, strip_d=FALSE),
                 c("IGHV1-69,IGHV1-69D", "IGHD1-1", "IGHJ1"))
    
    # Filtering non-localized genes
    nl_call <- c("IGHV3-NL1*01,IGHV3-30-3*01,IGHV3-30*01", 
                 "Homosap IGHV3-30*01 F,Homsap IGHV3-NL1*01 F",
                 "IGHV1-NL1*01")
    
    expect_equal(getAllele(nl_call, first=FALSE, omit_nl=TRUE),
                 c("IGHV3-30-3*01,IGHV3-30*01", "IGHV3-30*01", ""))
                 
    expect_equal(getGene(nl_call, first=FALSE, omit_nl=TRUE),
                 c("IGHV3-30-3,IGHV3-30", "IGHV3-30", ""))
    
    expect_equal(getFamily(nl_call, first=FALSE, omit_nl=TRUE),
                 c("IGHV3", "IGHV3", "IGHV1"))
    
    # Test for issue found in TIgGER
    # If there's no allele info,
    # getAllele should return the same input value
    # Before the change in allele_regex used by getAllele,
    # getAllele("TRAV38-2/DV8") returned TRAV38-2
    # which in TIgGER introduced NAs and an error, as TRAV38-2 was not a valid name
    # for the list of alleles, the valid name was TRAV38-2/DV8.
    # getAllele("IGHV1-01") returned IGHV1-01
    getAllele(getGene("IGHV1-01*02")) == getGene("IGHV1-01*02")
    getAllele(getGene("TRAV38-2/DV8*01")) == getGene("TRAV38-2/DV8*01")
})

### sortGenes ####

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

### groupGenes ####
test_that("groupGenes heavy only",{
    # make a data frame
    db <- data.frame(SEQUENCE_ID = c(1,2,3,4,5,6,7,8), 
                     V_CALL = c("IGHV1-1", "IGHV1-1,IGHV1-2", "IGHV1-2,IGHV1-3", "IGHV1-2", "IGHV1-3", "IGHV1-3", "IGHV1-1,IGHV1-2", "IGHV1-2"), 
                     J_CALL = c("IGHJ1", "IGHJ2", "IGHJ1", "IGHJ1", "IGHJ1", "IGHJ2", "IGHJ1", "IGHJ2"),
                     stringsAsFactors=F)
    # group VJ genes
    db <- groupGenes(db,
                     v_call="V_CALL",
                     j_call="J_CALL",
                     first=F)
    # test
    #expected <- c(2, 3, 2, 2, 2, 1, 2, 3)
    # same underlying spirit as before
    expected <- c("G1","G2","G1","G1","G1","G3","G1","G2")
    expect_equal(db$VJ_GROUP, expected)
})

test_that("groupGenes, single-cell mode, heavy and light", {
    
    # (in theory, should be 1 heavy per cell; but code-wise there is no restriction for groupGenes)
    load(file.path("..", "data-tests", "db_sc.rda")) # data_t1
    
    # View(data_t1[, c("V_CALL", "J_CALL", "LEN")])
    
    # manual deduction
    # 1 
    # IGHV3-11 IGHJ4 93 IGKV1-39 IGLJ6 30
    # 2 
    # IGHV3-7 IGHJ1 24 IGLV2-11 IGLJ6 33
    # IGHV3-6 IGHJ1 24 IGLV2-11 IGLJ6 33 (first=F)
    # IGHV1-4 IGHJ3 27 IGLV2-11 IGLJ6 33
    # IGHV1-4 IGHJ4 27 IGLV2-11 IGLJ6 33 (first=F)
    # IGHV3-7 IGHJ1 24 IGKV2-13 IGLJ3 57
    # IGHV3-6 IGHJ1 24 IGKV2-13 IGLJ3 57 (first=F)
    # IGHV1-4 IGHJ3 27 IGKV2-13 IGLJ3 57
    # IGHV1-4 IGHJ4 27 IGKV2-13 IGLJ3 57 (first=F)
    # 
    # 3 IGHV3-7 IGHJ5 90 IGLV3-30 IGLJ3 36
    # 4 IGHV3-7 IGHJ5 90 IGLV3-30 IGLJ3 36
    # 5 IGHV3-7 IGHJ5 90 IGLV2-20 IGLJ3 36
    #   IGHV3-7 IGHJ5 90 IGLV3-30 IGLJ3 36 (first=F)
    # 
    # 6 IGHV4-59 IGHJ1 84 IGKV1-27 IGKJ5 39
    #   IGHV4-59 IGHJ1 84 IGKV1-25 IGKJ3 60
    # 7 IGHV4-59 IGHJ1 84 IGKV1-27 IGKJ5 39
    # 8 IGHV4-59 IGHJ1 84 IGKV1-27 IGKJ5 39
    # 9 IGHV4-59 IGHJ1 84 IGKV1-27 IGKJ5 39
    # 
    # 10 IGHV4-59 IGHJ1 84 IGKV1-6 IGKJ4 42
    # 11 IGHV4-59 IGHJ1 84 IGKV1-6 IGKJ4 42
    
    # gg1 first=F
    # 1 by itself; 2 by itself; 3-5 together; 6-9 together; 10-11 together
    gg1_expect = c(rep("G2", 2), 
                   rep("G1", 3), 
                   rep("G3", 2+2+2), 
                   rep("G4", 3+2+2+2), 
                   rep("G5", 2+2))
    
    # gg2 first=T
    # 1 by itself; 2 by itself; 3-4 together; 5 by itself; 6-9 together; 10-11 together
    gg2_expect = c(rep("G2", 2),
                   rep("G1", 3), 
                   rep("G4", 2+2), 
                   rep("G3", 2), 
                   rep("G5", 3+2+2+2), 
                   rep("G6", 2+2))
    
    gg1 = groupGenes(data_t1, v_call="V_CALL", j_call="J_CALL", 
                     junc_len="LEN", cell_id="CELL_ID", locus="LOCUS", 
                     useOnlyIGH=FALSE, first=FALSE)
    
    gg2 = groupGenes(data_t1, v_call="V_CALL", j_call="J_CALL", 
                     junc_len="LEN", cell_id="CELL_ID", locus="LOCUS", 
                     useOnlyIGH=FALSE, first=TRUE)
                     
    expect_equal(gg1[["VJ_GROUP"]], gg1_expect)
    expect_equal(gg2[["VJ_GROUP"]], gg2_expect)
})

test_that("groupGenes, single-cell mode, heavy only", {
    
    # (in theory, should be 1 heavy per cell; but code-wise there is no restriction for groupGenes)
    load(file.path("..", "data-tests", "db_sc.rda")) # data_t1
    
    # manual deduction
    # 1 
    # IGHV3-11 IGHJ4 93 
    # 2 
    # IGHV3-7 IGHJ1 24 
    # IGHV3-6 IGHJ1 24  (first=F)
    # IGHV1-4 IGHJ3 27 
    # IGHV1-4 IGHJ4 27  (first=F)
    # 
    # 3 IGHV3-7 IGHJ5 90
    # 4 IGHV3-7 IGHJ5 90
    # 5 IGHV3-7 IGHJ5 90
    # 
    # 6 IGHV4-59 IGHJ1 84 
    # 7 IGHV4-59 IGHJ1 84 
    # 8 IGHV4-59 IGHJ1 84 
    # 9 IGHV4-59 IGHJ1 84 
    # 
    # 10 IGHV4-59 IGHJ1 84
    # 11 IGHV4-59 IGHJ1 84
    
    # 1 by itself, 2 by itself, 3-5 together, 6-11 together
    gg1_expect = c(rep("G2", 2), 
                   rep("G1", 3), 
                   rep("G3", 2+2+2), 
                   rep("G4", 3+2+2+2+2+2)) 
    
    gg1 = groupGenes(data_t1, v_call="V_CALL", j_call="J_CALL", 
                     junc_len="LEN", cell_id="CELL_ID", locus="LOCUS", 
                     useOnlyIGH=TRUE, first=FALSE)
    
    expect_equal(gg1[["VJ_GROUP"]], gg1_expect)
    
})
