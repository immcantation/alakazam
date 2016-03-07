file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
db <- readChangeoDb(file)

test_that("translateDNA", {
    expect_equal(translateDNA("ACTGACTCGA",trim=F),"TDS")
    expect_equal(translateDNA("ACTGACTCGA",trim=T),"D")
})

test_that("countOccurrences0",{
    expect_equal(countOccurrences("STSTSTS","STS"),2)
    #expect_equal(countOccurrences("STSTSTS","STS"),3)
})

test_that("gravy", {
    aa_seq <- "CARDRSTPWRRGIASTTVRTSW"
    expect_equal(gravy(aa_seq),-0.918, tolerance = .001)
})

test_that("bulk", {
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
    expect_equal(bulk(seq),c(14.46227,16.58857), tolerance = .001)
    
    data(aaindex, package="seqinr")
    x <- aaindex[["GRAR740103"]]$I
    # Rename the score vector to use single-letter codes
    names(x) <- translateStrings(names(x), ABBREV_AA)
    # Calculate average volume
    obs <- bulk(seq, bulkiness=x)
    expect_equal(obs, c(77.34091, 93.71429),tolerance=0.001)
})

test_that("polar", {
    # Default scale
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
    obs <- polar(seq)
    expect_equal(obs, c(8.55, 8.00),tolerance=0.001)
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
    obs <- aliphatic(seq)
    expect_equal(obs, c(0.4000000, NA, 0.4142857),tolerance=0.001)
    
})

test_that("charge", {
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
    # Normalized charge
    obs_norm <- charge(seq)
    expect_equal(obs_norm, c(0.1784859, 0.1425715),tolerance=0.001)
    
    # Unnormalized charge
    obs <- charge(seq, normalize=FALSE)
    expect_equal(obs, c(3.9266889, 0.9980008),tolerance=0.001)
    
    # Use the Murray et al, 2006 scores from the seqinr package
    data(pK, package="seqinr")
    x <- setNames(pK[["Murray"]], rownames(pK))
    # Calculate charge
    obs <- charge(seq, pK=x, normalize=FALSE)
    expect_equal(obs, c(3.8946562, 0.9977872),tolerance=0.001)
})

test_that("aminoAcidProperties", {
    seq_aa <- translateDNA(db$JUNCTION[1:5])
    
    junction_gravy <- gravy(seq_aa)
    junction_bulk <- bulk(seq_aa)
    junction_polar <- polar(seq_aa)
    junction_aliphatic <- aliphatic(seq_aa)
    junction_charge <- charge(seq_aa)
    
    junction_properties <- aminoAcidProperties(db[1:5,], seq="JUNCTION", nt=TRUE,
                                            trim=FALSE, label="JUNCTION")
    expect_equal(junction_gravy,junction_properties$JUNCTION_AA_GRAVY, tolerance = .001)
    expect_equal(junction_bulk,junction_properties$JUNCTION_AA_BULK, tolerance = .001)
    expect_equal(junction_polar,junction_properties$JUNCTION_AA_POLAR, tolerance = .001)
    expect_equal(junction_aliphatic,junction_properties$JUNCTION_AA_ALIPHATIC, tolerance = .001)
    expect_equal(junction_charge,junction_properties$JUNCTION_AA_CHARGE, tolerance = .001)
    
    
    data(aaindex, package="seqinr")
    h <- aaindex[["KIDA850101"]]$I
    # Rename the score vector to use single-letter codes
    names(h) <- translateStrings(names(h), ABBREV_AA)
    
    junction_gravy_h <- gravy(seq_aa, hydropathy =  h)
    junction_properties_h <- aminoAcidProperties(db[1:5,], seq="JUNCTION", nt=TRUE,
                                            trim=FALSE, label="JUNCTION",
                                            hydropathy = h)
    expect_equal(junction_gravy_h,junction_properties_h$JUNCTION_AA_GRAVY, tolerance = .001)
    
    junction_gravy_na <- gravy(c(NA,"NA","NULL"), hydropathy =  h)
    expect_equal(junction_gravy_na,c(NA,0.27,NA), tolerance = .001)
    
    db[1,"JUNCTION"] <- NA
    db[2,"JUNCTION"] <- "NA"
    db[3,"JUNCTION"] <- "NULL"
    db[4,"JUNCTION"] <- "NLL"
    junction_properties_na <- aminoAcidProperties(db[1:4,], seq="JUNCTION", nt=FALSE,
                                                 trim=FALSE, label="JUNCTION",
                                                 hydropathy = h)
    expect_equal(junction_properties_na$JUNCTION_AA_LENGTH,c(NA,2,NA,3))
    expect_equal(junction_properties_na$JUNCTION_AA_GRAVY, tolerance = .001,
                 c(NA,0.27,NA,-0.463))
    expect_equal(junction_properties_na$JUNCTION_AA_BULK,c(NA,12.16,NA,18.54),tolerance = .001)
    expect_equal(junction_properties_na$JUNCTION_AA_ALIPHATIC,c(NA,0.5,NA,2.6),tolerance = .001)
    expect_equal(junction_properties_na$JUNCTION_AA_POLARITY,c(NA,9.85,NA,7.13),tolerance = .001)
    expect_equal(junction_properties_na$JUNCTION_AA_CHARGE,c(NA,0,NA,0),tolerance = .001)
    expect_equal(junction_properties_na$JUNCTION_AA_BASIC,c(NA,0,NA,0),tolerance = .001)
    expect_equal(junction_properties_na$JUNCTION_AA_ACIDIC,c(NA,0,NA,0),tolerance = .001)
    expect_equal(isValidAASeq(db[1:4,"JUNCTION"]),c(F,T,F,T))
    expect_warning(aminoAcidProperties(db[1:4,], seq="JUNCTION", nt=FALSE,
                                       trim=FALSE, label="JUNCTION", property = "length"),
                   "2 sequences"
                   )
})

test_that("validate amino acid sequences" ,{
    seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVR--XX", "CARJ", "10")
    expect_equal(isValidAASeq(seq),c(T,T,F,F))

})

test_that("countPatterns", {
    seq <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
             "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
             "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
    patterns <- c("A", "V", "[LI]")
    names(patterns) <- c("ARG", "VAL", "ISO_LEU")
    obs <- countPatterns(seq, patterns, nt=TRUE, trim=TRUE, label="CDR3")
    expect_equal(obs$CDR3_ARG, c(0.1250000, 0.1111111, 0.1111111), tolerance=0.001)
    expect_equal(obs$CDR3_VAL, c(0, 0, 0))
    expect_equal(obs$CDR3_ISO_LEU, c(0, 0.1111111, 0), tolerance=0.001)
})