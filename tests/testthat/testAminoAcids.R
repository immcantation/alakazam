ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)

test_that("translateDNA", {
    expect_equal(translateDNA("ACTGACTCGA", trim=FALSE), "TDS")
    expect_equal(translateDNA("ACTGACTCGA", trim=TRUE), "D")
})

test_that("countOccurrences0",{
    expect_equal(countOccurrences("STSTSTS", "STS"), 2)
    #expect_equal(countOccurrences("STSTSTS","STS"), 3)
})

test_that("gravy", {
    aa_seq <- "CARDRSTPWRRGIASTTVRTSW"
    expect_equal(gravy(aa_seq), -0.918, tolerance=0.001)
})

test_that("bulk", {
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
    expect_equal(bulk(seq), c(14.46227, 16.58857), tolerance=0.001)
    
    data(aaindex, package="seqinr")
    x <- aaindex[["GRAR740103"]]$I
    # Rename the score vector to use single-letter codes
    names(x) <- translateStrings(names(x), ABBREV_AA)
    # Calculate average volume
    obs <- bulk(seq, bulkiness=x)
    expect_equal(obs, c(77.34091, 93.71429), tolerance=0.001)
})

test_that("polar", {
    # Default scale
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
    obs <- polar(seq)
    expect_equal(obs, c(8.55, 8.00), tolerance=0.001)
    # Use the Zimmerman et al, 1968 polarity scale from the seqinr package
    data(aaindex, package = "seqinr")
    x <- aaindex[["ZIMJ680103"]]$I
    # Rename the score vector to use single-letter codes
    names(x) <- translateStrings(names(x), ABBREV_AA)
    # Calculate polarity
    obs <- polar(seq, polarity=x)
    expect_equal(obs, c(14.94864,  8.86000),tolerance=0.001)
})


test_that("aliphatic", {
    seq <- c("CARDRSTPWRRGIASTTVRTSW", NA, "XXTQMYVRT")
    obs <- aliphatic(seq, normalize=TRUE)
    expect_equal(obs, c(0.4000000, NA, 0.4142857),tolerance=0.001)
    
})

test_that("charge", {
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
    # Unnormalized charge
    obs <- charge(seq, normalize=FALSE)
    expect_equal(obs, c(3.9266889, 0.9980008), tolerance=0.001)

    # Normalized charge
    obs_norm <- charge(seq, normalize=TRUE)
    expect_equal(obs_norm, c(0.1784859, 0.1425715),tolerance=0.001)
    
    # Use the Murray et al, 2006 scores from the seqinr package
    data(pK, package="seqinr")
    x <- setNames(pK[["Murray"]], rownames(pK))
    # Calculate charge
    obs <- charge(seq, pK=x, normalize=FALSE)
    expect_equal(obs, c(3.8946562, 0.9977872), tolerance=0.001)
})

test_that("aminoAcidProperties", {
    seq_aa <- translateDNA(db$JUNCTION[1:5])
    
    junction_gravy <- gravy(seq_aa)
    junction_bulk <- bulk(seq_aa)
    junction_polar <- polar(seq_aa)
    junction_aliphatic <- aliphatic(seq_aa)
    junction_charge <- charge(seq_aa)
    
    junction_properties <- aminoAcidProperties(db[1:5, ], seq="JUNCTION", 
                                               trim=FALSE, label="junction")
    expect_equal(junction_gravy, junction_properties$junction_aa_gravy, tolerance=0.001)
    expect_equal(junction_bulk, junction_properties$junction_aa_bulk, tolerance=0.001)
    expect_equal(junction_polar, junction_properties$junction_aa_polar, tolerance=0.001)
    expect_equal(junction_aliphatic, junction_properties$junction_aa_aliphatic, tolerance=0.001)
    expect_equal(junction_charge, junction_properties$junction_aa_charge, tolerance=0.001)
    
    
    data(aaindex, package="seqinr")
    h <- aaindex[["KIDA850101"]]$I
    # Rename the score vector to use single-letter codes
    names(h) <- translateStrings(names(h), ABBREV_AA)
    
    junction_gravy_h <- gravy(seq_aa, hydropathy =  h)
    junction_properties_h <- aminoAcidProperties(db[1:5,], seq="JUNCTION", 
                                            trim=FALSE, label="junction",
                                            hydropathy = h)
    expect_equal(junction_gravy_h, junction_properties_h$junction_aa_gravy, tolerance = .001)
    
    junction_gravy_na <- gravy(c(NA, "NA", "NULL"), hydropathy=h)
    expect_equal(junction_gravy_na, c(NA, 0.27, NA), tolerance=0.001)
    
    db[1,"JUNCTION"] <- NA
    db[2,"JUNCTION"] <- "NA"
    db[3,"JUNCTION"] <- "NULL"
    db[4,"JUNCTION"] <- "NLL"
    expect_warning(junction_properties_na <- aminoAcidProperties(db[1:4,], seq="JUNCTION", nt=FALSE,
                                                 trim=FALSE, label="junction",
                                                 hydropathy=h),
                   "Found 2 sequences with non valid amino acid symbols")
    expect_equal(junction_properties_na$junction_aa_length, c(NA, 2, NA, 3))
    expect_equal(junction_properties_na$junction_aa_gravy, tolerance=0.001,
                 c(NA,0.27,NA,-0.463))
    expect_equal(junction_properties_na$junction_aa_bulk,c(NA, 12.16, NA, 18.54), tolerance=0.001)
    expect_equal(junction_properties_na$junction_aa_aliphatic, c(NA, 0.5, NA, 2.6), tolerance=0.001)
    expect_equal(junction_properties_na$junction_aa_polarity, c(NA, 9.85, NA, 7.13), tolerance=0.001)
    expect_equal(junction_properties_na$junction_aa_charge, c(NA, 0, NA, 0), tolerance=0.001)
    expect_equal(junction_properties_na$junction_aa_basic, c(NA, 0, NA, 0), tolerance=0.001)
    expect_equal(junction_properties_na$junction_aa_acidic, c(NA, 0, NA, 0), tolerance=0.001)
    expect_equal(isValidAASeq(db$JUNCTION[1:4]), c(FALSE, TRUE, FALSE, TRUE))
    expect_warning(aminoAcidProperties(db[1:4, ], seq="JUNCTION", nt=FALSE,
                                       trim=FALSE, label="junction", property="length"),
                   "2 sequences")
    
    db <- data.frame(
        "junction"="ACTGACTCGA",
        "junction_aa"="TDS"
    )
    aa_prop_trim_F <- aminoAcidProperties(db, seq="junction_aa", nt=F, trim=F)
    expect_equal(aa_prop_trim_F$junction_aa_aa_length, 3)
    aa_prop_trim_T <- aminoAcidProperties(db, seq="junction_aa", nt=F, trim=T)
    expect_equal(aa_prop_trim_T$junction_aa_aa_length, 1)
})

test_that("validate amino acid sequences" ,{
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVR--XX", "CARJ", "10")
    expect_equal(isValidAASeq(seq), c(TRUE, TRUE, FALSE, FALSE))

})

test_that("countPatterns", {
    seq <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
             "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
             "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
    patterns <- c("A", "V", "[LI]")
    names(patterns) <- c("arg", "val", "iso_leu")
    obs <- countPatterns(seq, patterns, trim=TRUE, label="cdr3")
    expect_equal(obs$cdr3_arg, c(0.1250000, 0.1111111, 0.1111111), tolerance=0.001)
    expect_equal(obs$cdr3_val, c(0, 0, 0))
    expect_equal(obs$cdr3_iso_leu, c(0, 0.1111111, 0), tolerance=0.001)
})