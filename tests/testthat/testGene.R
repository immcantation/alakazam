ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)
expect_warning(
    db_gg <- readChangeoDb(file.path("..", "data-tests", "db_test.tsv")),
    "db_test.tsv is not in the Change-O format"
)

### countGenes ####

test_that("countGenes", {
    # Without copy numbers
    genes <- countGenes(db, gene = "V_CALL", groups = "SAMPLE", mode = "family")
    expect_equal(genes$seq_freq,
        c(0.76, 0.87, 0.24, 0.05, 0.05, 0.02, 0.01),
        tolerance = 0.001
    )

    genes <- countGenes(db, gene = "V_CALL", groups = "SAMPLE", mode = "gene")
    expect_equal(genes$seq_freq,
        c(
            0.41, 0.35, 0.24, 0.33, 0.26, 0.11, 0.10, 0.05, 0.05,
            0.04, 0.02, 0.02, 0.01, 0.01
        ),
        tolerance = 0.001
    )

    genes <- countGenes(db, gene = "V_CALL", groups = "SAMPLE", mode = "allele")
    expect_equal(genes$seq_freq,
        c(
            0.350, 0.325, 0.215, 0.330, 0.085, 0.150, 0.100, 0.100, 0.070,
            0.050, 0.050, 0.025, 0.040, 0.030, 0.020, 0.020,
            0.010, 0.010, 0.010, 0.010
        ),
        tolerance = 0.01
    )

    # With copy numbers and multiple groups
    genes <- countGenes(db,
        gene = "V_CALL", groups = c("SAMPLE", "ISOTYPE"),
        copy = "DUPCOUNT", mode = "family"
    )
    expect_equal(genes$seq_freq,
        c(
            0.78, 0.95, 0.95, 0.22, 0.63, 0.37, 0.60, 0.20, 0.40, 0.67,
            0.60, 0.80, 0.40, 0.33, 0.03, 0.05, 0.01, 0.01
        ),
        tolerance = 0.01
    )

    genes <- countGenes(db,
        gene = "V_CALL", groups = c("SAMPLE", "ISOTYPE"),
        copy = "DUPCOUNT", mode = "gene"
    )
    expect_equal(round(genes$seq_freq, 2)[1:12],
        c(
            0.61, 0.75, 0.22, 0.54, 0.40, 0.37, 0.17,
            0.20, 0.15, 0.20, 0.60, 0.20
        ),
        tolerance = 0.01
    )

    genes <- countGenes(db,
        gene = "V_CALL", groups = c("SAMPLE", "ISOTYPE"),
        copy = "DUPCOUNT", mode = "allele"
    )
    expect_equal(genes$seq_freq[1:12],
        c(0.61, 0.72, 0.16, 0.41, 0.40, 0.37, 0.06, 0.20, 0.08, 0.20, 0.50, 0.08),
        tolerance = 0.01
    )

    # Testing of fill

    # Without copy numbers
    genes <- countGenes(db, gene = "V_CALL", groups = "SAMPLE", mode = "family", fill = T)
    expect_equal(genes$seq_freq,
        c(0.02, 0.87, 0.05, 0.05, 0.01, 0, 0.76, 0, 0.24, 0),
        tolerance = 0.001
    )

    genes <- countGenes(db, gene = "V_CALL", groups = "SAMPLE", mode = "gene", fill = T)
    expect_equal(genes$seq_freq,
        c(
            0.02, 0.04, 0.01, 0.1, 0.11, 0.26, 0.33, 0.02, 0, 0.05, 0.05, 0.01, 0, 0,
            0, 0, 0.41, 0, 0, 0, 0.35, 0, 0.24, 0
        ),
        tolerance = 0.001
    )

    genes <- countGenes(db, gene = "V_CALL", groups = "SAMPLE", mode = "allele", fill = T)
    expect_equal(genes$seq_freq,
        c(
            0.02, 0.01, 0.03, 0.01, 0.1, 0.07, 0.04, 0.01, 0.1, 0.15, 0.33, 0.02, 0,
            0.05, 0.05, 0, 0.01, 0, 0, 0, 0, 0, 0.085, 0.325, 0, 0, 0, 0, 0, 0.35, 0, 0.215, 0.025, 0
        ),
        tolerance = 0.01
    )

    # With copy numbers and multiple groups
    genes <- countGenes(db,
        gene = "V_CALL", groups = c("SAMPLE", "ISOTYPE"),
        copy = "DUPCOUNT", mode = "family", fill = T
    )
    expect_equal(round(genes$seq_freq, 2)[1:12],
        c(0, 0.8, 0, 0.2, 0, 0, 0.6, 0.4, 0, 0, 0, 0.6),
        tolerance = 0.01
    )

    genes <- countGenes(db,
        gene = "V_CALL", groups = c("SAMPLE", "ISOTYPE"),
        copy = "DUPCOUNT", mode = "gene", fill = T
    )
    expect_equal(round(genes$seq_freq, 2)[1:12],
        c(0, 0, 0, 0.2, 0, 0, 0.6, 0, 0, 0, 0.2, 0),
        tolerance = 0.01
    )

    genes <- countGenes(db,
        gene = "V_CALL", groups = c("SAMPLE", "ISOTYPE"),
        copy = "DUPCOUNT", mode = "allele", fill = T
    )
    expect_equal(genes$seq_freq[1:12],
        c(0, 0, 0, 0, 0.2, 0, 0, 0, 0, 0, 0.6, 0),
        tolerance = 0.01
    )

    # Test how NAs are handled
    db_some_na <- data.frame(
        sequence_id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
        sample = c("S1", "S2", "S4", "S3", "S2", "S3", "S4", "S4", "S3", "S4"),
        v_call = c("IGHV1-1", NA, "IGHV1-1,IGHV1-2", "IGHV1-2,IGHV1-3", "IGHV1-2", "IGHV1-3", "IGHV1-3", "IGHV1-1,IGHV1-2", "IGHV1-2", NA),
        dupcount = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        clone = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4),
        stringsAsFactors = F
    )

    expect_warning(genes <- countGenes(db_some_na, gene = "v_call", mode = "gene"), NULL)
    expect_equal(genes$seq_freq,
        c(0.375, 0.375, 0.25),
        tolerance = 0.001
    )

    genes <- countGenes(db_some_na, gene = "v_call", mode = "gene", remove_na = FALSE)
    # A locus cannot be determined if the gene is NA. When determining the denominator (locus_count),
    # these sequences will be considered separately as their own "NA" locus.
    # Example:
    #  locus gene    seq_count locus_count seq_freq
    # <chr> <chr>       <int>       <int>    <dbl>
    # 1 IGH   IGHV1-1         3           8    0.375
    # 2 IGH   IGHV1-2         3           8    0.375
    # 3 IGH   IGHV1-3         2           8    0.25
    # 4 NA    NA              2           2    1
    expect_equal(genes$seq_freq,
        c(0.375, 0.375, 0.25, 1),
        tolerance = 0.001
    )

    genes <- countGenes(db_some_na, gene = "v_call", mode = "gene", copy = "clone", remove_na = FALSE)
    expect_equal(genes$copy_count,
        c(10, 7, 7, 6),
        tolerance = 0.001
    )

    genes <- countGenes(db_some_na, gene = "v_call", mode = "gene", groups = "sample", fill = TRUE, remove_na = FALSE)
    expect_equal(genes$seq_count,
        c(1, 0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 0, 1, 1),
        tolerance = 0.001
    )

    db_all_na <- data.frame(
        sequence_id = c(1, 2, 3),
        sample = c("S1", "S2", "S2"),
        v_call = c(NA, NA, NA),
        dupcount = c(1, 1, 1),
        clone = c(1, 1, 1),
        stringsAsFactors = F
    )

    expect_warning(genes <- countGenes(db_all_na, gene = "v_call", mode = "gene"), "The column v_call contains no data")
    expect_equal(nrow(genes),
        0,
        tolerance = 0.001
    )

    expect_warning(genes <- countGenes(db_all_na, gene = "v_call", mode = "gene", remove_na = FALSE), "The column v_call contains no data")
    expect_equal(genes$seq_count,
        3,
        tolerance = 0.001
    )

    expect_warning(
        genes <- countGenes(db_all_na, gene = "v_call", mode = "gene", groups = "sample", copy = "dupcount", fill = TRUE, remove_na = FALSE),
        "The column v_call contains no data"
    )
    expect_equal(genes$copy_count,
        c(1, 2),
        tolerance = 0.001
    )

    # testing addition of collapse and first flags

    genes <- countGenes(db, gene = "J_CALL", mode = "allele", first = FALSE)
    expect_equal(genes$seq_freq,
        c(0.283, 0.213, 0.183, 0.133, 0.123, 0.0233, 0.020, 0.010, 0.010),
        tolerance = 0.001
    )

    genes <- countGenes(db, gene = "V_CALL", mode = "gene", first = FALSE, collapse = FALSE)
    expect_equal(genes$seq_freq,
        c(
            0.31, 0.233, 0.17, 0.11, 0.0867, 0.0333, 0.0167, 0.0133, 0.00667,
            0.00667, 0.00667, 0.0033, 0.0033
        ),
        tolerance = 0.001
    )

    # test bulk percentages now being by locus
    # need data with mixed loci
    db_hl <- data.frame(
        sequence_id = c(1:30),
        sample = c(rep("S1", 15), rep("S2", 15)),
        v_call = c(
            "IGHV1-2,IGHV1-3", "IGHV1-2", "IGHV1-2", "IGHV1-3", "IGHV1-3", "IGHV1-3", "IGHV1-3", "IGHV2-1,IGHV2-2", "IGHV2-2", "IGHV2-2",
            "IGKV1-1", "IGKV1-2", "IGKV1-1,IGKV1-2", "IGLV1-1", "IGLV2-1",
            "IGHV1-2", "IGHV1-2", "IGHV3-1", "IGHV3-1", "IGHV2-1,IGHV2-2", "IGHV2-1", "IGHV2-2", "IGHV2-2", "IGHV2-2", "IGHV3-1",
            "IGKV1-1", "IGKV1-2", "IGKV2-1", "IGLV1-1", "IGLV1-1"
        ),
        j_call = c(
            "IGHJ1", "IGHJ1", "IGHJ1", "IGHJ1", "IGHJ2", "IGHJ2", "IGHJ2", "IGHJ2", "IGHJ2", "IGHJ2",
            "IGKJ1", "IGKJ1", "IGKJ1", "IGLJ1", "IGLJ1",
            "IGHJ1", "IGHJ1", "IGHJ1", "IGHJ1", "IGHJ2", "IGHJ2", "IGHJ2", "IGHJ2", "IGHJ2", "IGHJ2",
            "IGKJ1", "IGKJ2", "IGKJ2", "IGLJ1", "IGLJ1"
        ),
        dupcount = c(
            1, 1, 5, 1, 1, 1, 5, 1, 5, 1, 1, 4, 1, 4, 1,
            1, 1, 1, 1, 1, 2, 2, 4, 4, 1, 1, 4, 1, 4, 1
        ),
        clone = c(
            1, 1, 1, 1, 2, 2, 2, 3, 3, 3, NA, NA, NA, NA, NA,
            4, 4, 5, 5, 6, 6, 6, 6, 6, 7, NA, NA, NA, NA, NA
        ),
        stringsAsFactors = F
    )

    genes <- countGenes(db_hl, gene = "v_call", groups = "sample", mode = "family")
    expect_equal(genes$seq_freq,
        c(
            0.7, 0.5, 0.3, 1, 0.3, 0.2, 0.667, 1, 0.5, 0.5, 0.333
        ),
        tolerance = 0.001
    )

    genes <- countGenes(db_hl, gene = "v_call", groups = "sample", mode = "gene")
    expect_equal(genes$seq_freq,
        c(
            0.4, 0.3, 0.3, 0.3, 0.2, 0.667, 0.2, 0.2, 1, 0.1, 0.333, 0.5, 0.5, 0.333, 0.333, 0.333
        ),
        tolerance = 0.001
    )

    genes <- countGenes(db_hl, gene = "j_call", groups = "sample", mode = "gene")
    expect_equal(genes$seq_freq,
        c(
            0.6, 0.6, 0.4, 0.4, 1, 1, 0.667, 1, 0.333
        ),
        tolerance = 0.001
    )

    genes <- countGenes(db_hl, gene = "v_call", groups = "sample", mode = "gene", copy = "dupcount")
    expect_equal(genes$copy_freq,
        c(
            0.556, 0.364, 0.318, 0.273, 1, 0.667, 0.8, 0.667, 0.167, 0.167, 0.333, 0.111, 0.045, 0.2, 0.167, 0.167
        ),
        tolerance = 0.001
    )

    genes <- countGenes(db_hl, gene = "j_call", groups = "sample", mode = "gene", copy = "dupcount")
    expect_equal(genes$copy_freq,
        c(
            0.636, 0.778, 0.364, 1, 1, 0.833, 1, 0.222, 0.167
        ),
        tolerance = 0.001
    )

    # case with all clone IDs provided
    genes <- countGenes(db_hl %>% dplyr::filter(!is.na(clone)), gene = "v_call", groups = "sample", mode = "gene", clone = "clone")
    expect_equal(genes$clone_freq,
        c(
            0.5, 0.333, 0.333, 0.333, 0.25, 0.25
        ),
        tolerance = 0.001
    )

    # case with some clone IDs provided
    expect_warning(
        genes <- countGenes(db_hl, gene = "v_call", groups = "sample", mode = "gene", clone = "clone"),
        "Found 10 sequences without clonal assignments. These sequences will be removed before counting clone frequencies."
    )
    expect_equal(genes$clone_freq,
        c(
            0.5, 0.333, 0.333, 0.333, 0.25, 0.25
        ),
        tolerance = 0.001
    )

    expect_warning(
        genes <- countGenes(db_hl, gene = "j_call", groups = "sample", mode = "gene", clone = "clone"),
        "Found 10 sequences without clonal assignments. These sequences will be removed before counting clone frequencies."
    )
    expect_equal(genes$clone_freq,
        c(
            0.667, 0.5, 0.5, 0.333
        ),
        tolerance = 0.001
    )

    # test if all NA
    db_hl$clone_na <- rep(NA, nrow(db_hl))
    expect_error(
        genes <- countGenes(db_hl, gene = "v_call", groups = "sample", mode = "gene", clone = "clone_na"),
        "No clone IDs are present in the data."
    )
})

test_that("countGenes works for single-cell", {
    # Without copy numbers
    expect_warning(genes <- countGenes(db_gg, gene = "v_call", groups = "subject_id", mode = "family"), "Mixed bulk and single-cell data")
    expect_equal(round(genes$seq_freq, 2),
        c(0.83, 0.50, 0.11, 0.11, 0.06, 0.06),
        tolerance = 0.001
    )

    expect_warning(genes <- countGenes(db_gg, gene = "v_call", groups = "subject_id", mode = "gene"), "Mixed bulk and single-cell data")
    expect_equal(round(genes$seq_freq, 2),
        c(0.83, 0.50, 0.11, 0.11, 0.06, 0.06),
        tolerance = 0.001
    )

    added_cell <- data.frame(
        subject_id = c("S2", "S2", "S2"),
        v_call = c("IGHV1-1", "IGKV1-1", "IGKV1-1"),
        j_call = c("IGHJ1-1", "IGKJ1-1", "IGKJ1-1"),
        cell_id = c(30, 30, 30)
    )

    # Specific use case for identical light V gene within 1 cell that should be counted only once
    db_gg_multilight <- db_gg %>%
        dplyr::select(c("subject_id", "v_call", "j_call", "cell_id")) %>%
        dplyr::bind_rows(db_gg, added_cell)
    expect_warning(genes <- countGenes(db_gg_multilight, gene = "v_call", groups = "subject_id", mode = "gene"), "Mixed bulk and single-cell data")
    expect_equal(genes$seq_freq,
        c(0.8, 0.5, 0.1, 0.1, 0.05, 0.05, 1, 1),
        tolerance = 0.001
    )
})

### getSegment ####

test_that("getSegment", {
    # Light chain tests
    kappa_call <- c(
        "Homsap IGKV1D-39*01 F,Homsap IGKV1-39*02 F,Homsap IGKV1-39*01",
        "Homsap IGKJ5*01 F"
    )

    expect_equal(
        getAllele(kappa_call),
        c("IGKV1-39*01", "IGKJ5*01")
    )
    expect_equal(
        getAllele(kappa_call, first = FALSE),
        c("IGKV1-39*01,IGKV1-39*02", "IGKJ5*01")
    )
    expect_equal(
        getAllele(kappa_call, first = FALSE, strip_d = FALSE),
        c("IGKV1D-39*01,IGKV1-39*02,IGKV1-39*01", "IGKJ5*01")
    )

    expect_equal(
        getGene(kappa_call),
        c("IGKV1-39", "IGKJ5")
    )
    expect_equal(
        getGene(kappa_call, first = FALSE),
        c("IGKV1-39", "IGKJ5")
    )
    expect_equal(
        getGene(kappa_call, first = FALSE, strip_d = FALSE),
        c("IGKV1D-39,IGKV1-39", "IGKJ5")
    )
    expect_equal(
        getGene(kappa_call, first = FALSE, strip_d = TRUE),
        c("IGKV1-39", "IGKJ5")
    )

    expect_equal(
        getFamily(kappa_call),
        c("IGKV1", "IGKJ5")
    )
    expect_equal(
        getFamily(kappa_call, first = FALSE),
        c("IGKV1", "IGKJ5")
    )
    expect_equal(
        getFamily(kappa_call, first = FALSE, collapse = FALSE),
        c("IGKV1,IGKV1,IGKV1", "IGKJ5")
    )
    expect_equal(
        getFamily(kappa_call, first = FALSE, strip_d = FALSE),
        c("IGKV1D,IGKV1", "IGKJ5")
    )

    expect_equal(
        getLocus(kappa_call, first = FALSE, strip_d = FALSE, collapse = TRUE),
        c("IGK", "IGK")
    )
    expect_equal(
        getLocus(kappa_call, first = FALSE, strip_d = FALSE, collapse = FALSE),
        c("IGK,IGK,IGK", "IGK")
    )

    expect_equal(
        getChain(kappa_call, first = FALSE, strip_d = FALSE, collapse = TRUE),
        c("VL", "VL")
    )
    expect_equal(
        getChain(kappa_call, first = FALSE, strip_d = FALSE, collapse = FALSE),
        c("VL,VL,VL", "VL")
    )

    # Heavy chain tests
    heavy_call <- c(
        "Homsap IGHV1-69*01 F,Homsap IGHV1-69D*01 F",
        "Homsap IGHD1-1*01 F",
        "Homsap IGHJ1*01 F"
    )

    expect_equal(
        getAllele(heavy_call, first = FALSE),
        c("IGHV1-69*01", "IGHD1-1*01", "IGHJ1*01")
    )
    expect_equal(
        getAllele(heavy_call, first = FALSE, strip_d = FALSE),
        c("IGHV1-69*01,IGHV1-69D*01", "IGHD1-1*01", "IGHJ1*01")
    )

    expect_equal(
        getGene(heavy_call, first = FALSE),
        c("IGHV1-69", "IGHD1-1", "IGHJ1")
    )
    expect_equal(
        getGene(heavy_call, first = FALSE, strip_d = FALSE),
        c("IGHV1-69,IGHV1-69D", "IGHD1-1", "IGHJ1")
    )
    expect_equal(
        getGene(heavy_call, first = FALSE, strip_d = TRUE),
        c("IGHV1-69", "IGHD1-1", "IGHJ1")
    )

    expect_equal(
        getLocus(heavy_call, first = FALSE, strip_d = FALSE, collapse = TRUE),
        c("IGH", "IGH", "IGH")
    )
    expect_equal(
        getLocus(heavy_call, first = FALSE, strip_d = FALSE, collapse = FALSE),
        c("IGH,IGH", "IGH", "IGH")
    )

    expect_equal(
        getChain(heavy_call, first = FALSE, strip_d = FALSE, collapse = TRUE),
        c("VH", "VH", "VH")
    )
    expect_equal(
        getChain(heavy_call, first = FALSE, strip_d = FALSE, collapse = FALSE),
        c("VH,VH", "VH", "VH")
    )

    # Filtering non-localized genes
    nl_call <- c(
        "IGHV3-NL1*01,IGHV3-30-3*01,IGHV3-30*01",
        "Homosap IGHV2-30*01 F,Homsap IGHV3-NL1*01 F",
        "IGHV1-NL1*01"
    )

    expect_equal(
        getAllele(nl_call, first = FALSE, omit_nl = TRUE),
        c("IGHV3-30-3*01,IGHV3-30*01", "IGHV2-30*01", "")
    )

    expect_equal(
        getGene(nl_call, first = FALSE, omit_nl = TRUE),
        c("IGHV3-30-3,IGHV3-30", "IGHV2-30", "")
    )

    expect_equal(
        getFamily(nl_call, first = FALSE, omit_nl = TRUE),
        c("IGHV3", "IGHV2", "")
    )

    expect_equal(
        getLocus(nl_call, first = FALSE, omit_nl = TRUE),
        c("IGH", "IGH", "")
    )

    expect_equal(
        getChain(nl_call, first = FALSE, omit_nl = TRUE),
        c("VH", "VH", "")
    )

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
    gene <- countGenes(data.frame("c_call" = "IGHD", stringsAsFactors = F), gene = "c_call", mode = "gene")
    asis <- countGenes(data.frame("c_call" = "IGHD", stringsAsFactors = F), gene = "c_call", mode = "asis")
    expect_equal(gene[["gene"]], "IGHD")
    expect_equal(asis[["gene"]], "IGHD")

    # Test for IMGT temporary designation nomenclature
    tmp_call <- c("IGHV9S3*01", "IGKV10S12*01")
    expect_equal(
        getAllele(tmp_call),
        c("IGHV9S3*01", "IGKV10S12*01")
    )
    expect_equal(
        getGene(tmp_call),
        c("IGHV9S3", "IGKV10S12")
    )
    expect_equal(
        getFamily(tmp_call),
        c("IGHV9", "IGKV10")
    )


    # Test Mus musculus
    # ADEGMC
    mm_alleles <- c(
        "IGHA*01",
        "IGHD*01", "IGHD5-1-1*01",
        "IGHE*01",
        "IGHG1*01", "IGHG2A*01",
        "IGHM*01",
        "IGKC*01"
    )
    mm_genes <- c(
        "IGHA",
        "IGHD", "IGHD5-1-1",
        "IGHE",
        "IGHG1", "IGHG2A",
        "IGHM",
        "IGKC"
    )
    mm_families <- c(
        "IGHA",
        "IGHD", "IGHD5",
        "IGHE",
        "IGHG1", "IGHG2A",
        "IGHM",
        "IGKC"
    )
    expect_equal(getGene(mm_alleles), mm_genes)
    expect_equal(getFamily(mm_alleles), mm_families)


    # Test for TCR
    expect_equal(
        getGene(c("TRBV29-1*01", "TRAV12-1*01")),
        c("TRBV29-1", "TRAV12-1")
    )
    expect_equal(
        getFamily(c("TRBV29-1*01", "TRAV12-1*01")),
        c("TRBV29", "TRAV12")
    )
})

### sortGenes ####

test_that("sortGenes", {
    genes <- c(
        "IGHV1-69D*01", "IGHV1-69*01", "IGHV4-38-2*01", "IGHV1-69-2*01",
        "IGHV2-5*01", "IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05",
        "IGHV1-2", "IGHV1-2*02", "IGHV1-69*02"
    )
    # Sort genes by name
    expected <- c(
        "IGHV1-2", "IGHV1-2*01,IGHV1-2*05", "IGHV1-2*02",
        "IGHV1-69*01", "IGHV1-69D*01", "IGHV1-69*02", "IGHV1-69-2*01",
        "IGHV1-NL1*01", "IGHV2-5*01", "IGHV4-38-2*01"
    )
    sorted <- sortGenes(genes)
    expect_equal(sorted, expected)

    # Sort genes by position in the locus
    expected_locus <- c(
        "IGHV1-NL1*01", "IGHV1-69-2*01", "IGHV1-69*01",
        "IGHV1-69D*01", "IGHV1-69*02", "IGHV4-38-2*01",
        "IGHV2-5*01", "IGHV1-2", "IGHV1-2*01,IGHV1-2*05",
        "IGHV1-2*02"
    )
    sorted_locus <- sortGenes(genes, method = "pos")
    expect_equal(sorted_locus, expected_locus)
})

### groupGenes ####

test_that("groupGenes heavy only", {
    # make a data frame
    db <- data.frame(
        SEQUENCE_ID = c(1, 2, 3, 4, 5, 6, 7, 8),
        V_CALL = c("IGHV1-1", "IGHV1-1,IGHV1-2", "IGHV1-2,IGHV1-3", "IGHV1-2", "IGHV1-3", "IGHV1-3", "IGHV1-1,IGHV1-2", "IGHV1-2"),
        J_CALL = c("IGHJ1", "IGHJ2", "IGHJ1", "IGHJ1", "IGHJ1", "IGHJ2", "IGHJ1", "IGHJ2"),
        stringsAsFactors = F
    )
    db$LOCUS <- getLocus(db$V_CALL)
    # group VJ genes
    db <- groupGenes(db,
        v_call = "V_CALL",
        j_call = "J_CALL",
        locus = "LOCUS",
        first = F
    )
    # test
    # expected <- c(2, 3, 2, 2, 2, 1, 2, 3)
    # same underlying spirit as before
    expected <- c("G1", "G2", "G1", "G1", "G1", "G3", "G1", "G2")
    expect_equal(db$vj_group, expected)
})


test_that("groupGenes, when 1 row ", {
    # Tests fix for rowSums when 1 row
    db <- structure(list(
        sample_id = "UW0343__S", sequence_id = "assemble12_0",
        subject_id = "S1", v_call = "IGHV3-30*02,IGHV3-30-5*02",
        j_call = "IGHJ6*03", junction = "TGTGCGAAAGTCCCCGTGGGGACTGCCTCTTACATGGACGCCTGG",
        cell_id = NA_character_, single_cell = FALSE, locus = "IGH",
        junction_length = 45L
    ), row.names = c(NA, -1L), class = c("tbl_df", "tbl", "data.frame"))
    ft <- groupGenes(db,
        junc_len = "junction_length",
        cell_id = "cell_id", locus = "locus", only_heavy = T,
        first = TRUE
    )
    ff <- groupGenes(db,
        junc_len = "junction_length",
        cell_id = "cell_id", locus = "locus", only_heavy = T,
        first = FALSE
    )
    expect_equal(ft, ff)
})


test_that("groupGenes all cell_id values are NA", {
    db <- data.frame(
        subject_id = c("S1", "S1", "S1", "S1", "S1", "S1", "S1", "S1"),
        v_call = c("IGHV1-1*01", "IGHV1-1*01", "IGHV1-2*01", "IGHV1-1*01,IGHV1-2*01", "IGHV1-2*01", "IGKV1-1*01", "IGKV1-1*01","IGKV1-1*01"),
        j_call = c("IGHJ2*01", "IGHJ1*01", "IGHJ1*01", "IGHJ1*01", "IGHJ1*01", "IGKJ1*01", "IGKJ1*01", "IGKJ1*01"),
        junction = c("TGTAAAAAATGG", "TGTAAAAAATGG", "TGTAAAACCTGG", "TGTAAACCCTGG", "TGTAAACCCTGG", "TGTCCCCCCTGG", "TGTCCCCCCTGG", "TGTCCCCCCTGG"),
        locus = c("IGH", "IGH", "IGH", "IGH", "IGH", "IGK", "IGK", "IGK"),
        cell_id = c(NA, NA, NA, NA, NA, NA, NA, NA),
        junction_length = 12
    )

    expect_warning(
        groupGenes(db,
                        v_call = "v_call", j_call = "j_call", junc_len = NULL,
                        locus = "locus", only_heavy = TRUE,
                        first = FALSE),
        "A cell_id column was found in the data, but was not specified. All values are NA."
    )    
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
        singleCellValidation(data_t1, cell_id = "CELL_ID", locus = "LOCUS"),
        regexp = "Only one heavy chain is allowed per cell"
    )
    expect_error(
        gg1 <- groupGenes(data_t1,
            v_call = "V_CALL", j_call = "J_CALL",
            junc_len = "LEN", cell_id = "CELL_ID", locus = "LOCUS",
            only_heavy = TRUE, first = FALSE
        ),
        regexp = "Only one heavy chain is allowed per cell"
    )

    # 1 by itself, 2 by itself, 3-5 together, 6-11 together
    gg1_expect <- c(
        rep("G1", 2),
        # rep("G1", 4), # Cell B removed
        rep("G2", 2 + 2 + 2),
        rep("G3", 2 + 3 + 2 + 2 + 2 + 2)
    )

    gg1 <- groupGenes(data_t1 %>% dplyr::filter(CELL_ID != "B"),
        v_call = "V_CALL", j_call = "J_CALL",
        junc_len = "LEN", cell_id = "CELL_ID", locus = "LOCUS",
        only_heavy = TRUE, first = FALSE
    )

    expect_equal(gg1[["vj_group"]], gg1_expect)
})

test_that("groupGenes, mixed bulk and single cell", {
    # TODO: add other sc cell with different L to test only_heavy T/F and fix tests
    db <- data.frame(
        subject_id = c("S1", "S1", "S1", "S1", "S1", "S1", "S1", "S1"),
        v_call = c("IGHV1-1*01", "IGHV1-1*01", "IGHV1-2*01", "IGHV1-1*01,IGHV1-2*01", "IGHV1-2*01", "IGKV1-1*01", "IGKV1-1*01","IGKV1-1*01"),
        j_call = c("IGHJ2*01", "IGHJ1*01", "IGHJ1*01", "IGHJ1*01", "IGHJ1*01", "IGKJ1*01", "IGKJ1*01", "IGKJ1*01"),
        junction = c("TGTAAAAAATGG", "TGTAAAAAATGG", "TGTAAAACCTGG", "TGTAAACCCTGG", "TGTAAACCCTGG", "TGTCCCCCCTGG", "TGTCCCCCCTGG", "TGTCCCCCCTGG"),
        locus = c("IGH", "IGH", "IGH", "IGH", "IGH", "IGK", "IGK", "IGK"),
        cell_id = c(1, 2, 3, NA, NA, 1, NA, 4),
        junction_length = 12
    )

    # subset to single-cell sequences
    # expect warning, one cell has only light chain and will be assigned group NA
    expect_warning(
        d <- groupGenes(db %>% dplyr::filter(!is.na(cell_id)),
                        v_call = "v_call", j_call = "j_call", junc_len = NULL,
                        cell_id = "cell_id", locus = "locus", only_heavy = TRUE,
                        first = FALSE),
                        pattern="NA(s) found in one or more of "
        )
    expect_equal(d[["vj_group"]], c("G2", "G1", "G3", "G2", NA))

    expect_warning(
        d_oh <- groupGenes(db %>% dplyr::filter(!is.na(cell_id)),
            v_call = "v_call", j_call = "j_call", junc_len = NULL,
            cell_id = "cell_id", locus = "locus", only_heavy = FALSE,
            first = FALSE
        ),
        "only_heavy = FALSE is deprecated"
    )

    expect_warning(
        d_sl <- groupGenes(db %>% dplyr::filter(!is.na(cell_id)),
            v_call = "v_call", j_call = "j_call", junc_len = NULL,
            cell_id = "cell_id", locus = "locus", only_heavy = TRUE,
            first = FALSE, split_light = TRUE
        ),
        "split_light = TRUE is deprecated"
    )
    expect_equal(d[["vj_group"]], d_oh[["vj_group"]], d_sl[["vj_group"]])

    # bulk data
    d_bulk <- groupGenes(db %>% dplyr::filter(is.na(cell_id)),
        v_call = "v_call", j_call = "j_call", junc_len = NULL,
        cell_id = "cell_id", locus = "locus", only_heavy = TRUE,
        first = FALSE
    )

    # mixed data
    expect_warning(
        d_mixed <- groupGenes(db,
                              v_call = "v_call", j_call = "j_call", junc_len = NULL,
                              cell_id = "cell_id", locus = "locus", only_heavy = TRUE,
                              first = FALSE
        ),
        pattern="NA(s) found in one or more of "
    )
    expect_equal(d_mixed[["vj_group"]], c("G2", "G1", "G1", "G1", "G1", "G2", NA, NA))
    # TODO SSNN 7/16/2024
    # Run tests with the same data used in scoper

    # Bulk mode cell_id=NULL
    ## first=T
    gg <- groupGenes(db_gg %>% select(-cell_id), cell_id = NULL, first = T)
    expect_equal(
        gg[["vj_group"]],
        gg[["expected_group_cell_id-null_first-T"]]
    )
    ## first=F
    gg <- groupGenes(db_gg %>% select(-cell_id), cell_id = NULL, first = F)
    expect_equal(
        gg[["vj_group"]],
        gg[["expected_group_cell_id-null_first-F"]]
    )
})


#### AIRR-format migration tests ####

test_that("countGenes, AIRR-format migration", {
    db_c <- readChangeoDb(file.path("..", "data-tests", "ExampleDb.gz"))
    expect_warning(
        db_a <- readChangeoDb(
            file.path("..", "data-tests", "ExampleDb_airr.gz")
        ),
        regexp = "airr::read_rearrangement"
    )

    # 1
    genes_c <- countGenes(db_c, gene = "V_CALL", groups = "SAMPLE", mode = "allele")
    genes_a <- countGenes(db_a, gene = "v_call", groups = "sample", mode = "allele")

    expect_true(all(genes_c == genes_a))

    # 2
    genes_c <- countGenes(db_c,
        gene = "V_CALL", groups = c("SAMPLE", "ISOTYPE"),
        copy = "DUPCOUNT", mode = "family"
    )
    genes_a <- countGenes(db_a,
        gene = "v_call", groups = c("sample", "isotype"),
        copy = "duplicate_count", mode = "family"
    )

    expect_true(all(genes_c == genes_a))

    # 3
    genes_c <- countGenes(db_c,
        gene = "V_CALL", groups = "SAMPLE",
        mode = "allele", fill = TRUE
    )
    genes_a <- countGenes(db_a,
        gene = "v_call", groups = "sample",
        mode = "allele", fill = TRUE
    )

    # genes_c uses "SAMPLE", genes_a uses "sample"
    # make colnames same case for comparison
    colnames(genes_c) <- tolower(colnames(genes_c))
    expect_true(all.equal(genes_c, genes_a[rownames(genes_c), names(genes_c)]))
})

test_that("groupGenes, AIRR-format migration", {
    # 1
    db_c <- readChangeoDb(file.path("..", "data-tests", "ExampleDb.gz"))
    expect_warning(
        db_a <- readChangeoDb(file.path("..", "data-tests", "ExampleDb_airr.gz")),
        regexp = "airr::read_rearrangement"
    )

    db_c$LOCUS <- getLocus(db_c$V_CALL)

    newDb_c <- groupGenes(db_c, v_call = "V_CALL", j_call = "J_CALL", locus = "LOCUS")
    newDb_a <- groupGenes(db_a, v_call = "v_call", j_call = "j_call")

    # newDb_c and newDb_a not directly comparable b/c diff ncol
    expect_true(all(newDb_c[["vj_group"]] == newDb_a[["vj_group"]]))

    rm(db_c, db_a, newDb_c, newDb_a)

    # 2
    # Cell with CELL_ID==B has 2 heavy and 2 light chains
    load(file.path("..", "data-tests", "db_sc.rda")) # data_t1, data_t2
    load(file.path("..", "data-tests", "db_sc_airr.rda")) # data_t1_airr, data_t2_airr

    newDb_c <- groupGenes(data_t1 %>% dplyr::filter(CELL_ID != "B"),
        v_call = "V_CALL", j_call = "J_CALL",
        junc_len = "LEN", cell_id = "CELL_ID", locus = "LOCUS",
        only_heavy = TRUE, first = FALSE
    )
    newDb_a <- groupGenes(data_t1_airr %>% dplyr::filter(cell_id != "B"),
        v_call = "v_call", j_call = "j_call",
        junc_len = "len", cell_id = "cell_id", locus = "locus",
        only_heavy = TRUE, first = FALSE
    )

    expect_true(all(newDb_c[["vj_group"]] == newDb_a[["vj_group"]]))
})

test_that("groupGenes, bulk mode, NA heavy chain V call is removed with a warning", {
    # Regression test: an NA in v_call made the one-to-one annotation check
    # evaluate if (NA) and stop with "missing value where TRUE/FALSE needed"
    # before the NA rows could be removed.
    db_na <- data.frame(
        sequence_id = c("s1", "s2", "s3"),
        v_call = c("IGHV1-2*01", NA, "IGHV1-2*01"),
        j_call = c("IGHJ4*02", "IGHJ4*02", "IGHJ4*02"),
        junction_length = c(42, 42, 42),
        locus = "IGH",
        stringsAsFactors = FALSE
    )
    expect_warning(
        res <- groupGenes(db_na, v_call = "v_call", j_call = "j_call"),
        "NA\\(s\\) found in one or more of \\{ v_call, j_call \\} columns. 1 sequence\\(s\\) removed"
    )
    expect_equal(res$sequence_id, c("s1", "s3"))
    expect_equal(length(unique(res$vj_group)), 1)
    expect_false(any(is.na(res$vj_group)))
})

#### groupLightGenes ####

# Build a single-cell data set from one row per chain:
#   cell_id, locus, v_call, j_call, sequence_alignment, clone_id
makeLightDb <- function(...) {
    rows <- list(...)
    db <- data.frame(
        sequence_id = paste0("s", seq_along(rows)),
        cell_id = vapply(rows, `[[`, character(1), 1),
        locus = vapply(rows, `[[`, character(1), 2),
        v_call = vapply(rows, `[[`, character(1), 3),
        j_call = vapply(rows, `[[`, character(1), 4),
        sequence_alignment = vapply(rows, `[[`, character(1), 5),
        clone_id = vapply(rows, `[[`, character(1), 6),
        stringsAsFactors = FALSE
    )
    db$v_call[db$v_call == "NA"] <- NA
    db$sequence_alignment[db$sequence_alignment == "NA"] <- NA
    db$junction_length <- ifelse(db$locus == "IGH", 42, 30)
    db
}
KV1 <- "IGKV1-39*01"; KJ1 <- "IGKJ1*01"
KV3 <- "IGKV3-20*01"; KJ2 <- "IGKJ2*01"
LV2 <- "IGLV2-14*01"; LJ2 <- "IGLJ2*01"
HV <- "IGHV1-2*01"; HJ <- "IGHJ4*02"

# group of the chain(s) of `cell` with `locus` in result `res`
grp <- function(res, cell, locus = "IGH") {
    unique(res$vj_group_light[res$cell_id == cell & res$locus == locus])
}

test_that("groupLightGenes, paired and unpaired heavy chains", {
    db <- makeLightDb(
        # clone 1: two light chain groups (IGK and IGL)
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"),
        c("c1", "IGK", KV1, KJ1, "TTT", "1"),
        c("c2", "IGH", HV, HJ, "AAAAAAAAAC", "1"),
        c("c2", "IGK", KV1, KJ1, "TTT", "1"),
        c("c3", "IGH", HV, HJ, "CCCCCCCCCC", "1"),
        c("c3", "IGL", LV2, LJ2, "TTT", "1"),
        # unpaired heavy chains, closest to c1/c2 and to c3
        c("c4", "IGH", HV, HJ, "AAAAAAAAAT", "1"),
        c("c5", "IGH", HV, HJ, "CCCCCCCCCA", "1"),
        # clone 2: no light chains at all
        c("c6", "IGH", HV, HJ, "AAAAAAAAAA", "2"),
        c("c7", "IGH", HV, HJ, "AAAAAAAAAA", "2")
    )
    res <- groupLightGenes(db)

    expect_equal(nrow(res), nrow(db))
    expect_equal(res$sequence_id, db$sequence_id)
    # heavy chains inherit the group of their cell's light chain
    for (cell in c("c1", "c2", "c3")) {
        expect_equal(grp(res, cell, "IGH"), grp(res, cell, ifelse(cell == "c3", "IGL", "IGK")))
    }
    expect_equal(grp(res, "c1"), grp(res, "c2"))
    expect_false(grp(res, "c1") == grp(res, "c3"))
    # unpaired heavy chains go to the group of their nearest paired heavy chain
    expect_equal(grp(res, "c4"), grp(res, "c1"))
    expect_equal(grp(res, "c5"), grp(res, "c3"))
    # clones without light chains stay NA
    expect_true(all(is.na(res$vj_group_light[res$clone_id == "2"])))
    # labels have the groupGenes format
    expect_true(all(grepl("^G[0-9]+$", res$vj_group_light[res$clone_id == "1"])))
})

test_that("groupLightGenes, cells with multiple light chains link groups", {
    db <- makeLightDb(
        # c1 has an IGK and an IGL light chain, linking the IGK V1 and IGL groups
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"),
        c("c1", "IGK", KV1, KJ1, "TTT", "1"), 
        c("c1", "IGL", LV2, LJ2, "TTT", "1"),
        c("c2", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c2", "IGK", KV1, KJ1, "TTT", "1"),
        c("c3", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c3", "IGL", LV2, LJ2, "TTT", "1"),
        # different light chain gene, not linked
        c("c4", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c4", "IGK", KV3, KJ2, "TTT", "1")
    )
    res <- groupLightGenes(db)

    linked <- unique(res$vj_group_light[res$cell_id %in% c("c1", "c2", "c3")])
    expect_length(linked, 1)
    expect_false(grp(res, "c4") == linked)
    expect_equal(grp(res, "c4", "IGK"), grp(res, "c4"))
    # merged groups are renumbered without gaps
    expect_setequal(unique(res$vj_group_light), c("G1", "G2"))
})

test_that("groupLightGenes, ties between unpaired heavy chain neighbours", {
    # Tied neighbours in different groups: the larger group wins.
    # u1 is at distance 5 from c1 (group of c1/c2, size 2) and from c3 (size 1)
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "TTT", "1"),
        c("c2", "IGH", HV, HJ, "GGGGGGGGGG", "1"), 
        c("c2", "IGK", KV1, KJ1, "TTT", "1"),
        c("c3", "IGH", HV, HJ, "CCCCCCCCCC", "1"), 
        c("c3", "IGL", LV2, LJ2, "TTT", "1"),
        c("u1", "IGH", HV, HJ, "AAAAACCCCC", "1")
    )
    res <- groupLightGenes(db)
    expect_equal(grp(res, "u1"), grp(res, "c1"))

    # Tied neighbours in equally sized groups: the lowest numbered group wins.
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "TTT", "1"),
        c("c3", "IGH", HV, HJ, "CCCCCCCCCC", "1"), 
        c("c3", "IGL", LV2, LJ2, "TTT", "1"),
        c("u1", "IGH", HV, HJ, "AAAAACCCCC", "1")
    )
    res <- groupLightGenes(db)
    candidates <- c(grp(res, "c1"), grp(res, "c3"))
    lowest <- candidates[which.min(as.integer(sub("^G", "", candidates)))]
    expect_equal(grp(res, "u1"), lowest)
})

test_that("groupLightGenes, light chains with missing gene calls", {
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "TTT", "1"),
        c("c3", "IGH", HV, HJ, "CCCCCCCCCC", "1"), 
        c("c3", "IGL", LV2, LJ2, "TTT", "1"),
        # c2's only light chain has no V call: c2 is treated as unpaired
        c("c2", "IGH", HV, HJ, "AAAAAAAAAC", "1"), 
        c("c2", "IGK", "NA", KJ1, "TTT", "1")
    )
    expect_warning(
        res <- groupLightGenes(db),
        "1 light chain sequence\\(s\\) with NA or empty values"
    )
    # the heavy chain is placed by distance, and the light chain inherits its group
    expect_equal(grp(res, "c2"), grp(res, "c1"))
    expect_equal(grp(res, "c2", "IGK"), grp(res, "c1"))
})

test_that("groupLightGenes, NA heavy chain sequence_alignment is an error", {
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "TTT", "1"),
        c("c2", "IGH", HV, HJ, "NA", "1")
    )
    expect_error(
        groupLightGenes(db),
        "1 heavy chain sequence\\(s\\) have NA or empty values in the sequence_alignment column"
    )
    db$sequence_alignment[3] <- ""
    expect_error(groupLightGenes(db), "NA or empty values in the sequence_alignment column")

    # missing light chain alignments are not used, so they are allowed
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "NA", "1")
    )
    expect_silent(res <- groupLightGenes(db))
    expect_equal(grp(res, "c1"), "G1")
})

test_that("groupLightGenes, input validation", {
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "TTT", "1")
    )
    # cell_id = NULL: warning and all NA
    expect_warning(
        res <- groupLightGenes(db, cell_id = NULL),
        "cell_id not specified, but 'cell_id' found in data"
    )
    expect_true(all(is.na(res$vj_group_light)))
    expect_warning(
        res <- groupLightGenes(db[, names(db) != "cell_id"], cell_id = NULL),
        "cell_id not specified\\. Heavy chain"
    )
    expect_true(all(is.na(res$vj_group_light)))
    # cell_id specified but column missing: error
    expect_error(groupLightGenes(db[, names(db) != "cell_id"]),
                 "column is not found in the data")

    # missing required column
    expect_error(groupLightGenes(db[, names(db) != "sequence_alignment"]),
                 "One or more required columns were not found")

    # invalid locus
    db_bad <- db
    db_bad$locus[2] <- "XYZ"
    expect_error(groupLightGenes(db_bad), "invalid loci annotations")

    # nproc must be a positive integer
    expect_error(groupLightGenes(db, nproc = 0), "nproc must be a positive integer")
    expect_error(groupLightGenes(db, nproc = 1.5), "nproc must be a positive integer")
    expect_error(groupLightGenes(db, nproc = "2"), "nproc must be a positive integer")
    expect_error(groupLightGenes(db, nproc = c(1, 2)), "nproc must be a positive integer")

    # factors are converted to character
    db_f <- db
    db_f$v_call <- factor(db_f$v_call)
    db_f$locus <- factor(db_f$locus)
    res <- groupLightGenes(db_f)
    expect_type(res$v_call, "character")
    expect_type(res$locus, "character")
    expect_equal(res$vj_group_light, c("G1", "G1"))
})

test_that("groupLightGenes, nproc > 1 gives the same result as nproc = 1", {
    skip_on_os("windows")
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "TTT", "1"),
        c("c2", "IGH", HV, HJ, "CCCCCCCCCC", "1"), 
        c("c2", "IGL", LV2, LJ2, "TTT", "1"),
        c("c3", "IGH", HV, HJ, "CCCCCCCCCA", "1"),
        c("c4", "IGH", HV, HJ, "AAAAAAAAAA", "2"), 
        c("c4", "IGK", KV3, KJ2, "TTT", "2"),
        c("c5", "IGH", HV, HJ, "AAAAAAAAAA", "3"), 
        c("c5", "IGK", KV1, KJ1, "TTT", "3"),
        c("c6", "IGH", HV, HJ, "AAAAAAAAAA", "4")
    )
    expect_equal(groupLightGenes(db, nproc = 2), groupLightGenes(db, nproc = 1))
})

test_that("groupLightGenes, reports warnings and errors from every clone", {
    skip_on_os("windows")
    db <- makeLightDb(
        c("c1", "IGH", HV, HJ, "AAAAAAAAAA", "1"), 
        c("c1", "IGK", KV1, KJ1, "TTT", "1"),
        c("c2", "IGH", HV, HJ, "AAAAAAAAAA", "2"), 
        c("c2", "IGK", KV1, KJ1, "TTT", "2"),
        # clone 3 has no light chains, so groupGenes is not called for it
        c("c3", "IGH", HV, HJ, "AAAAAAAAAA", "3")
    )
    real_groupGenes <- groupGenes

    # warnings raised inside a clone reach the caller, also from forked
    # workers, once per message
    local_mocked_bindings(groupGenes = function(...) {
        warning("mock groupGenes warning")
        real_groupGenes(...)
    })
    for (n in c(1, 2)) {
        warns <- character(0)
        res <- withCallingHandlers(
            groupLightGenes(db, nproc = n),
            warning = function(w) {
                warns <<- c(warns, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        expect_equal(warns, "mock groupGenes warning [clone(s): 1, 2]", info = paste("nproc =", n))
        expect_equal(res$vj_group_light, c("G1", "G1", "G1", "G1", NA), info = paste("nproc =", n))
    }

    # an error in any clone stops the run and names the clone
    local_mocked_bindings(groupGenes = function(data, ...) {
        if (any(data$v_call == KV1)) stop("mock groupGenes error")
        real_groupGenes(data, ...)
    })
    for (n in c(1, 2)) {
        expect_error(
            groupLightGenes(db, nproc = n),
            "groupLightGenes failed for clone\\(s\\):\n  1: mock groupGenes error\n  2: mock groupGenes error",
            info = paste("nproc =", n)
        )
    }
})