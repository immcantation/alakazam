file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
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

test_that("getProperties", {
    seq_aa <- db$JUNCTION[1:5]
    junction_gravy <- gravy(seq_aa, nt=TRUE)
    junction_properties <- regionProperties(db[1:5,], seq="JUNCTION", nt=TRUE,
                                            trim=FALSE, label="JUNCTION")
    expect_equal(junction_gravy,junction_properties$JUNCTION_AA_GRAVY, tolerance = .001)
    
    data(aaindex)
    h <- aaindex[["KIDA850101"]]$I
    # Rename the score vector to use single-letter codes
    names(h) <- translateStrings(names(h), AA_TRANS)
    
    seq_aa <- db$JUNCTION[1:5]
    junction_gravy <- gravy(seq_aa, hydropathy =  h, nt=TRUE)
    junction_properties <- regionProperties(db[1:5,], seq="JUNCTION", nt=TRUE,
                                            trim=FALSE, label="JUNCTION",
                                            hydropathy = h)
    expect_equal(junction_gravy,junction_properties$JUNCTION_AA_GRAVY, tolerance = .001)
    
})