ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)
expect_warning(db_gg <- readChangeoDb(file.path("..", "data-tests", "db_test.tsv")),
               "db_test.tsv is not in the Change-O format")

### countGenes ####

test_that("countGenes", {
    # Without copy numbers
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="family")
    expect_equal(genes$seq_freq, 
                 c(0.76, 0.87, 0.24, 0.05, 0.05, 0.02, 0.01),
                 tolerance=0.001)
    
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="gene")
    expect_equal(genes$seq_freq, 
                 c(0.41, 0.35, 0.24, 0.33, 0.26, 0.11, 0.10, 0.05, 0.05,
                   0.04, 0.02, 0.02, 0.01, 0.01),
                 tolerance=0.001)
    
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="allele")
    expect_equal(genes$seq_freq, 
                 c(0.350, 0.325, 0.215, 0.330, 0.085, 0.150, 0.100, 0.100, 0.070,
                   0.050, 0.050, 0.025, 0.040, 0.030, 0.020, 0.020,
                   0.010, 0.010, 0.010, 0.010),
                 tolerance=0.01)
    
    # With copy numbers and multiple groups
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="family")
    expect_equal(genes$seq_freq, 
                 c(0.78, 0.95, 0.95, 0.22, 0.63, 0.37, 0.60, 0.20, 0.40, 0.67,
                   0.60, 0.80, 0.40, 0.33, 0.03, 0.05, 0.01, 0.01),
                 tolerance=0.01)
    
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="gene")
    expect_equal(round(genes$seq_freq,2)[1:12], 
                 c(0.61, 0.75, 0.22, 0.54, 0.40, 0.37, 0.17,
                   0.20, 0.15, 0.20, 0.60, 0.20),
                 tolerance=0.01)
    
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="allele")
    expect_equal(genes$seq_freq[1:12], 
                 c(0.61, 0.72, 0.16, 0.41, 0.40, 0.37, 0.06, 0.20, 0.08, 0.20, 0.50, 0.08),
                 tolerance=0.01)

    # Testing of fill

    # Without copy numbers
    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="family", fill = T)
    expect_equal(genes$seq_freq, 
                 c(0.02,0.87,0.05,0.05,0.01,0,0.76,0,0.24,0),
                 tolerance=0.001)

    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="gene", fill = T)
    expect_equal(genes$seq_freq, 
                 c(0.02,0.04,0.01,0.1,0.11,0.26,0.33,0.02,0,0.05,0.05,0.01,0,0,
                    0,0,0.41,0,0,0,0.35,0,0.24,0),
                 tolerance=0.001)

    genes <- countGenes(db, gene="V_CALL", groups="SAMPLE", mode="allele", fill = T)
    expect_equal(genes$seq_freq, 
                 c(0.02,0.01,0.03,0.01,0.1,0.07,0.04,0.01,0.1,0.15,0.33,0.02,0,
                    0.05,0.05,0,0.01,0,0,0,0,0,0.085,0.325,0,0,0,0,0,0.35,0,0.215,0.025,0),
                 tolerance=0.01)

    # With copy numbers and multiple groups
    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="family", fill = T)
    expect_equal(round(genes$seq_freq,2)[1:12], 
                 c(0,0.8,0,0.2,0,0,0.6,0.4,0,0,0,0.6),
                 tolerance=0.01)

    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="gene", fill = T)
    expect_equal(round(genes$seq_freq,2)[1:12], 
                 c(0,0,0,0.2,0,0,0.6,0,0,0,0.2,0),
                 tolerance=0.01)

    genes <- countGenes(db, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                        copy="DUPCOUNT", mode="allele", fill = T)
    expect_equal(genes$seq_freq[1:12], 
                 c(0,0,0,0,0.2,0,0,0,0,0,0.6,0),
                 tolerance=0.01)

    # Test how NAs are handled
    db_some_na <- data.frame(sequence_id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
                             sample = c("S1", "S2", "S4", "S3", "S2", "S3", "S4", "S4", "S3", "S4"),
                             v_call = c("IGHV1-1", NA, "IGHV1-1,IGHV1-2", "IGHV1-2,IGHV1-3", "IGHV1-2", "IGHV1-3", "IGHV1-3", "IGHV1-1,IGHV1-2", "IGHV1-2", NA), 
                             dupcount = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
                             clone = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4), 
                             stringsAsFactors=F)
    
    expect_warning(genes <- countGenes(db_some_na, gene="v_call", mode="gene"), NULL)
    expect_equal(genes$seq_freq, 
                 c(0.375, 0.375, 0.25),
                 tolerance=0.001)
    
    genes <- countGenes(db_some_na, gene="v_call", mode="gene", remove_na = FALSE)
    expect_equal(genes$seq_freq, 
                 c(0.3, 0.3, 0.2, 0.2),
                 tolerance=0.001)
    
    genes <- countGenes(db_some_na, gene="v_call", mode="gene", copy = "clone", remove_na = FALSE)
    expect_equal(genes$copy_count, 
                 c(10, 7, 7, 6),
                 tolerance=0.001)
    
    genes <- countGenes(db_some_na, gene="v_call", mode="gene", groups = "sample", fill = TRUE, remove_na = FALSE)
    expect_equal(genes$seq_count, 
                 c(1, 0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 0, 1, 1),
                 tolerance=0.001)
    
    db_all_na <- data.frame(sequence_id = c(1, 2, 3),
                            sample = c("S1", "S2", "S2"),
                            v_call = c(NA, NA, NA), 
                            dupcount = c(1, 1, 1), 
                            clone = c(1, 1, 1), 
                            stringsAsFactors=F)
    
    expect_warning(genes <- countGenes(db_all_na, gene="v_call", mode="gene"), "The column v_call contains no data")
    expect_equal(nrow(genes), 
                 0,
                 tolerance=0.001)
    
    expect_warning(genes <- countGenes(db_all_na, gene="v_call", mode="gene", remove_na=FALSE), "The column v_call contains no data")
    expect_equal(genes$seq_count, 
                 3,
                 tolerance=0.001)
    
    expect_warning(genes <- countGenes(db_all_na, gene="v_call", mode="gene", groups = "sample", copy = "dupcount", fill = TRUE,  remove_na = FALSE), 
                   "The column v_call contains no data")
    expect_equal(genes$copy_count, 
                 c(1, 2),
                 tolerance=0.001)
    
})

test_that("countGenes works for single-cell", {
    # Without copy numbers
    expect_warning(genes <- countGenes(db_gg, gene="v_call", groups="subject_id", mode="family"), "Mixed bulk and single-cell data")
    expect_equal(round(genes$seq_freq,2), 
                 c(0.83, 0.50, 0.11, 0.11, 0.06, 0.06),
                 tolerance=0.001)
    
    expect_warning(genes <- countGenes(db_gg, gene="v_call", groups="subject_id", mode="gene"), "Mixed bulk and single-cell data")
    expect_equal(round(genes$seq_freq,2), 
                 c(0.83, 0.50, 0.11, 0.11, 0.06, 0.06),
                 tolerance=0.001)
    
    added_cell <- data.frame(
        subject_id = c("S2","S2","S2"),
        v_call = c("IGHV1-1","IGKV1-1","IGKV1-1"),
        j_call = c("IGHJ1-1","IGKJ1-1","IGKJ1-1"),
        cell_id = c(30,30,30)
    )

    # Specific use case for identical light V gene within 1 cell that should be counted only once
    db_gg_multilight <- db_gg %>%
                        dplyr::select(c("subject_id","v_call","j_call","cell_id")) %>%
                        dplyr::bind_rows(db_gg, added_cell)
    expect_warning(genes <- countGenes(db_gg_multilight, gene="v_call", groups="subject_id", mode="gene"), "Mixed bulk and single-cell data")
    expect_equal(genes$seq_freq, 
                 c(0.8, 0.5, 0.1, 0.1, 0.05, 0.05, 1,1),
                 tolerance=0.001)
})

### getSegment ####

test_that("getSegment", {
    # Light chain tests
    kappa_call <- c("Homsap IGKV1D-39*01 F,Homsap IGKV1-39*02 F,Homsap IGKV1-39*01",
                    "Homsap IGKJ5*01 F")
    
    expect_equal(getAllele(kappa_call), 
                 c("IGKV1-39*01", "IGKJ5*01"))
    expect_equal(getAllele(kappa_call, first=FALSE),
                 c("IGKV1-39*01,IGKV1-39*02", "IGKJ5*01"))
    expect_equal(getAllele(kappa_call, first=FALSE, strip_d=FALSE),
                 c("IGKV1D-39*01,IGKV1-39*02,IGKV1-39*01", "IGKJ5*01" ))
    
    expect_equal(getGene(kappa_call), 
                 c("IGKV1-39", "IGKJ5" ))
    expect_equal(getGene(kappa_call, first=FALSE),
                 c("IGKV1-39", "IGKJ5"))
    expect_equal(getGene(kappa_call, first=FALSE, strip_d=FALSE),
                 c("IGKV1D-39,IGKV1-39", "IGKJ5"))
    expect_equal(getGene(kappa_call, first=FALSE, strip_d=TRUE),
                 c("IGKV1-39", "IGKJ5"))    
    
    expect_equal(getFamily(kappa_call), 
                 c("IGKV1", "IGKJ5"))
    expect_equal(getFamily(kappa_call, first=FALSE), 
                 c("IGKV1", "IGKJ5"))
    expect_equal(getFamily(kappa_call, first=FALSE, collapse=FALSE),
                 c("IGKV1,IGKV1,IGKV1", "IGKJ5"))
    expect_equal(getFamily(kappa_call, first=FALSE, strip_d=FALSE),
                 c("IGKV1D,IGKV1", "IGKJ5"))

    expect_equal(getLocus(kappa_call, first=FALSE, strip_d=FALSE, collapse=TRUE),
                 c("IGK", "IGK"))    
    expect_equal(getLocus(kappa_call, first=FALSE, strip_d=FALSE, collapse=FALSE),
                 c("IGK,IGK,IGK", "IGK"))
    
    expect_equal(getChain(kappa_call, first=FALSE, strip_d=FALSE, collapse=TRUE),
                 c("VL", "VL"))
    expect_equal(getChain(kappa_call, first=FALSE, strip_d=FALSE, collapse=FALSE),
                 c("VL,VL,VL", "VL"))
    
    # Heavy chain tests
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
    expect_equal(getGene(heavy_call, first=FALSE, strip_d=TRUE),
                 c("IGHV1-69", "IGHD1-1", "IGHJ1"))    
    
    expect_equal(getLocus(heavy_call, first=FALSE, strip_d=FALSE, collapse=TRUE),
                 c("IGH", "IGH", "IGH"))    
    expect_equal(getLocus(heavy_call, first=FALSE, strip_d=FALSE, collapse=FALSE),
                 c("IGH,IGH", "IGH", "IGH"))    

    expect_equal(getChain(heavy_call, first=FALSE, strip_d=FALSE, collapse=TRUE),
                 c("VH", "VH", "VH"))        
    expect_equal(getChain(heavy_call, first=FALSE, strip_d=FALSE, collapse=FALSE),
                 c("VH,VH", "VH", "VH"))    
    
    # Filtering non-localized genes
    nl_call <- c("IGHV3-NL1*01,IGHV3-30-3*01,IGHV3-30*01", 
                 "Homosap IGHV2-30*01 F,Homsap IGHV3-NL1*01 F",
                 "IGHV1-NL1*01")
    
    expect_equal(getAllele(nl_call, first=FALSE, omit_nl=TRUE),
                 c("IGHV3-30-3*01,IGHV3-30*01", "IGHV2-30*01", ""))
                 
    expect_equal(getGene(nl_call, first=FALSE, omit_nl=TRUE),
                 c("IGHV3-30-3,IGHV3-30", "IGHV2-30", ""))
    
    expect_equal(getFamily(nl_call, first=FALSE, omit_nl=TRUE),
                 c("IGHV3", "IGHV2", ""))
    
    expect_equal(getLocus(nl_call, first=FALSE, omit_nl=TRUE),
                 c("IGH", "IGH", ""))
    
    expect_equal(getChain(nl_call, first=FALSE, omit_nl=TRUE),
                 c("VH", "VH", ""))
    
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
    
    # Test that getSegment works well with c_call and strip_d=TRUE
    # and it doesn't convert IGHD into IGH
    expect_equal(getGene("IGHD"), "IGHD")
    expect_equal(getAllele("IGHD"), "IGHD")
    gene <- countGenes(data.frame("c_call"="IGHD", stringsAsFactors = F),gene="c_call", mode="gene")
    asis <- countGenes(data.frame("c_call"="IGHD", stringsAsFactors = F),gene="c_call", mode="asis")
    expect_equal(gene[['gene']], "IGHD")
    expect_equal(asis[['gene']], "IGHD")
    
    # Test for IMGT temporary designation nomenclature
    tmp_call <- c("IGHV9S3*01", "IGKV10S12*01")
    expect_equal(getAllele(tmp_call),
                 c("IGHV9S3*01", "IGKV10S12*01"))
    expect_equal(getGene(tmp_call),
                 c("IGHV9S3", "IGKV10S12"))
    expect_equal(getFamily(tmp_call),
                 c("IGHV9", "IGKV10"))
    
    
    # Test Mus musculus
    # ADEGMC
    mm_alleles <- c("IGHA*01", 
                    "IGHD*01","IGHD5-1-1*01",
                    "IGHE*01",
                    "IGHG1*01","IGHG2A*01",
                    "IGHM*01",
                    "IGKC*01")
    mm_genes <- c("IGHA",
                  "IGHD","IGHD5-1-1",
                  "IGHE",
                  "IGHG1","IGHG2A",
                  "IGHM",
                  "IGKC")
    mm_families <- c("IGHA",
                  "IGHD","IGHD5",
                  "IGHE",
                  "IGHG1","IGHG2A",
                  "IGHM",
                  "IGKC")    
    expect_equal(getGene(mm_alleles),mm_genes)
    expect_equal(getFamily(mm_alleles),mm_families)
    
    
    # Test for TCR
    expect_equal(getGene(c("TRBV29-1*01","TRAV12-1*01")),
                 c("TRBV29-1","TRAV12-1"))
    expect_equal(getFamily(c("TRBV29-1*01","TRAV12-1*01")),
                 c("TRBV29","TRAV12"))    
                 
})

### sortGenes ####

test_that("sortGenes", {
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

test_that("groupGenes heavy only", {
    # make a data frame
    db <- data.frame(SEQUENCE_ID = c(1,2,3,4,5,6,7,8), 
                     V_CALL = c("IGHV1-1", "IGHV1-1,IGHV1-2", "IGHV1-2,IGHV1-3", "IGHV1-2", "IGHV1-3", "IGHV1-3", "IGHV1-1,IGHV1-2", "IGHV1-2"), 
                     J_CALL = c("IGHJ1", "IGHJ2", "IGHJ1", "IGHJ1", "IGHJ1", "IGHJ2", "IGHJ1", "IGHJ2"),
                     stringsAsFactors=F)
    db$LOCUS <- getLocus(db$V_CALL)
    # group VJ genes
    db <- groupGenes(db,
                     v_call="V_CALL",
                     j_call="J_CALL",
                     locus="LOCUS",
                     first=F)
    # test
    #expected <- c(2, 3, 2, 2, 2, 1, 2, 3)
    # same underlying spirit as before
    expected <- c("G1","G2","G1","G1","G1","G3","G1","G2")
    expect_equal(db$vj_group, expected)

})


test_that("groupGenes, when 1 row ", {
    # Tests fix for rowSums when 1 row
    db <- structure(list(sample_id = "UW0343__S", sequence_id = "assemble12_0", 
                   subject_id = "S1", v_call = "IGHV3-30*02,IGHV3-30-5*02", 
                   j_call = "IGHJ6*03", junction = "TGTGCGAAAGTCCCCGTGGGGACTGCCTCTTACATGGACGCCTGG", 
                   cell_id = NA_character_, single_cell = FALSE, locus = "IGH", 
                   junction_length = 45L), row.names = c(NA, -1L
                   ), class = c("tbl_df", "tbl", "data.frame"))
    ft <- groupGenes(db, junc_len="junction_length",
                  cell_id="cell_id", locus="locus", only_heavy=T,
                  first=TRUE)
    ff <- groupGenes(db, junc_len="junction_length",
                     cell_id="cell_id", locus="locus", only_heavy=T,
                     first=FALSE)
    expect_equal(ft,ff)
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
    
    # Cell with CELL_ID==B has 2 heavy and 2 light chains
    # expect error
    expect_error(
        singleCellValidation(data_t1, cell_id="CELL_ID",locus="LOCUS"),
        regexp = "Only one heavy chain is allowed per cell"
    )
    expect_error(
        gg1 <- groupGenes(data_t1, v_call="V_CALL", j_call="J_CALL", 
                          junc_len="LEN", cell_id="CELL_ID", locus="LOCUS", 
                          only_heavy=TRUE, first=FALSE),
        regexp = "Only one heavy chain is allowed per cell"
    )
    
    # 1 by itself, 2 by itself, 3-5 together, 6-11 together
    gg1_expect = c(rep("G1", 2), 
                   # rep("G1", 4), # Cell B removed
                   rep("G2", 2+2+2), 
                   rep("G3", 2+3+2+2+2+2)) 
    
    gg1 = groupGenes(data_t1 %>% dplyr::filter(CELL_ID != "B"), v_call="V_CALL", j_call="J_CALL", 
                     junc_len="LEN", cell_id="CELL_ID", locus="LOCUS", 
                     only_heavy=TRUE, first=FALSE)
    
    expect_equal(gg1[["vj_group"]], gg1_expect)
    
})

test_that("groupGenes, mixed bulk and single cell", {

    # TODO: add other sc cell with different L to test only_heavy T/F and fix tests
    db <- data.frame(
        subject_id=c("S1","S1","S1","S1","S1", "S1", "S1"),
        v_call=c("IGHV1-1*01","IGHV1-1*01","IGHV1-2*01","IGHV1-1*01,IGHV1-2*01","IGHV1-2*01", "IGKV1-1*01", "IGKV1-1*01"),
        j_call=c("IGHJ2*01","IGHJ1*01","IGHJ1*01","IGHJ1*01","IGHJ1*01","IGKJ1*01", "IGKJ1*01"),
        junction=c("TGTAAAAAATGG","TGTAAAAAATGG","TGTAAAACCTGG","TGTAAACCCTGG","TGTAAACCCTGG","TGTCCCCCCTGG","TGTCCCCCCTGG"),
        locus=c("IGH","IGH","IGH","IGH","IGH","IGK","IGK"),
        cell_id=c(1,2,3,NA,NA,1,NA),
        junction_length=12
    )


    # cell_id='cell_id', only_heavy=T
    expect_warning(
        c <- groupGenes(db, v_call="v_call", j_call="j_call", junc_len=NULL,
               cell_id="cell_id", locus="locus", only_heavy=TRUE,
               first=FALSE),
        "Mixed single cell and bulk data was indicated"
    )
    
    # TODO
    # Should match a, because mix bulk-sc -> bulk mode
    # expect_equal(c,a)

    # subset to single-cell sequences
    d <- groupGenes(db %>% dplyr::filter(!is.na(cell_id)), v_call="v_call", j_call="j_call", junc_len=NULL,
                    cell_id="cell_id", locus="locus", only_heavy=TRUE,
                    first=FALSE)
    expect_equal(d[['vj_group']], c("G2","G1","G3","G2"))

    # TODO SSNN 7/16/2024
    # Run tests with the same data used in scoper

    # Bulk mode cell_id=NULL
    ## first=T
    gg <- groupGenes(db_gg %>% select(-cell_id), cell_id=NULL, first=T)
    expect_equal(gg[["vj_group"]],
                 gg[["expected_group_cell_id-null_first-T"]])
    ## first=F
    gg <- groupGenes(db_gg %>% select(-cell_id), cell_id=NULL, first=F)
    expect_equal(gg[["vj_group"]], 
                 gg[["expected_group_cell_id-null_first-F"]])
    
    # Single cell mode cell_id="cell_id" mixed with bulk
    ## first=TRUE
    
    ### only_heavy=TRUE
    # TODO fix this. groupGenes now removes bulk sequences
    expect_warning(
        gg <- groupGenes(db_gg, cell_id="cell_id", only_heavy=TRUE, first=TRUE),
        "Mixed single cell and bulk data was indicated"
    )
    # gg[["vj_group"]]
    # gg[["expected_group_cell_id_first-T_only_heavy-T"]]
    
       
    ## first=FALSE
    #TODO
    
    ### only_heavy=TRUE
    #TODO
    

})


#### AIRR-format migration tests ####

test_that("countGenes, AIRR-format migration", {
    
    db_c <- readChangeoDb(file.path("..", "data-tests", "ExampleDb.gz"))
    expect_warning(
        db_a <- readChangeoDb(
            file.path("..", "data-tests", "ExampleDb_airr.gz")),
        regexp = "airr::read_rearrangement"
    )
    
    # 1
    genes_c <- countGenes(db_c, gene="V_CALL", groups="SAMPLE", mode="allele")
    genes_a <- countGenes(db_a, gene="v_call", groups="sample", mode="allele")
    
    expect_true(all(genes_c==genes_a))
    
    # 2
    genes_c <- countGenes(db_c, gene="V_CALL", groups=c("SAMPLE", "ISOTYPE"), 
                          copy="DUPCOUNT", mode="family")
    genes_a <- countGenes(db_a, gene="v_call", groups=c("sample", "isotype"), 
                          copy="duplicate_count", mode="family")
    
    expect_true(all(genes_c==genes_a))
    
    # 3
    genes_c <- countGenes(db_c, gene="V_CALL", groups="SAMPLE", 
                          mode="allele", fill=TRUE)
    genes_a <- countGenes(db_a, gene="v_call", groups="sample", 
                          mode="allele", fill=TRUE)
    
    expect_true(all(genes_c==genes_a))
    
})

test_that("groupGenes, AIRR-format migration", {

    # 1 
    db_c <- readChangeoDb(file.path("..", "data-tests", "ExampleDb.gz"))
    expect_warning(
        db_a <- readChangeoDb(file.path("..", "data-tests", "ExampleDb_airr.gz")),
        regexp="airr::read_rearrangement"
    )
    
    db_c$LOCUS <- getLocus(db$V_CALL)

    newDb_c <- groupGenes(db_c, v_call="V_CALL", j_call="J_CALL", locus="LOCUS")
    newDb_a <- groupGenes(db_a, v_call="v_call", j_call="j_call")
    
    # newDb_c and newDb_a not directly comparable b/c diff ncol
    expect_true(all(newDb_c[["vj_group"]]==newDb_a[["vj_group"]]))
    
    rm(db_c, db_a, newDb_c, newDb_a)

    # 2 
    # Cell with CELL_ID==B has 2 heavy and 2 light chains
    load(file.path("..", "data-tests", "db_sc.rda")) # data_t1, data_t2
    load(file.path("..", "data-tests", "db_sc_airr.rda")) # data_t1_airr, data_t2_airr

    newDb_c <- groupGenes(data_t1 %>% dplyr::filter(CELL_ID != "B"), 
                          v_call="V_CALL", j_call="J_CALL", 
                          junc_len="LEN", cell_id="CELL_ID", locus="LOCUS", 
                          only_heavy=TRUE, first=FALSE)
    newDb_a <- groupGenes(data_t1_airr %>% dplyr::filter(cell_id != "B"),
                          v_call="v_call", j_call="j_call", 
                          junc_len="len", cell_id="cell_id", locus="locus", 
                          only_heavy=TRUE, first=FALSE)
    
    expect_true(all(newDb_c[["vj_group"]]==newDb_a[["vj_group"]]))
    
})
