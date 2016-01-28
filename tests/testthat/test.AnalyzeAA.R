file <- system.file("extdata", "changeo_demo.gz", package="alakazam")
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

test_that("aminoAcidProperties", {
    seq_aa <- translateDNA(db$JUNCTION[1:5])
    
    junction_gravy <- gravy(seq_aa)
    junction_properties <- aminoAcidProperties(db[1:5,], seq="JUNCTION", nt=TRUE,
                                            trim=FALSE, label="JUNCTION")
    expect_equal(junction_gravy,junction_properties$JUNCTION_AA_GRAVY, tolerance = .001)
    
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
                                       trim=FALSE, label="JUNCTION",
                                       hydropathy = h, property = "length"),
                   "2 sequences"
                   )
})

