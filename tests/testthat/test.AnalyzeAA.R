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
