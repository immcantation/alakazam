file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

test_that("seqDist",{
    # Ungapped examples
    expect_equal(seqDist("ATGGC", "ATGGG"), 1)
    expect_equal(seqDist("ATGGC", "ATG??"), 2)
    
    # Gaps will be treated as Ns with a gap=0 distance matrix
    expect_equal(
        seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=0)),
        0)
    
    # Gaps will be treated as universally non-matching characters with gap=1
    expect_equal(
        seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=1)),
        2)
    
    # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
    expect_equal(
        seqDist("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1)),
        1)
    
    # Gaps of equivalent run lengths are not counted as gaps
    expect_equal(
        seqDist("ATG-C", "ATG-C", dist_mat=getDNAMatrix(gap=-1)),
        0)
    
    # Overlapping runs of gap characters are counted as a single gap
    expect_equal(
        seqDist("ATG-C", "AT--C", dist_mat=getDNAMatrix(gap=-1)),
        1) 
    
    expect_equal(
        seqDist("A-GGC", "AT--C", dist_mat=getDNAMatrix(gap=-1)),
        1)
    
    expect_equal(
        seqDist("AT--C", "AT--C", dist_mat=getDNAMatrix(gap=-1)),
        0)
    
    # Discontiguous runs of gap characters each count as separate gaps
    expect_equal(
        seqDist("-TGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1)),
        2)
})

test_that("pairwiseDist", {
    # Gaps will be treated as Ns with a gap=0 distance matrix
    obs <- pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
                 dist_mat=getDNAMatrix(gap=0))
    expect_equal(obs,
                 matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0),ncol=4),
                 check.attributes=F)
    # Gaps will be treated as universally non-matching characters with gap=1
    obs <- pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
                 dist_mat=getDNAMatrix(gap=1))
    expect_equal(obs,
                 matrix(c(0, 1, 1, 2, 1, 0, 0, 3, 1, 0, 0, 3, 2, 3, 3, 0),ncol=4),
                 check.attributes=F)
    
    # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
    obs <- pairwiseDist(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
                 dist_mat=getDNAMatrix(gap=-1))
    expect_equal(obs,
                 matrix(c(0, 1, 1, 1, 1, 0, 0, 2, 1, 0, 0, 2, 1, 2, 2, 0),ncol=4),
                 check.attributes=F)
})

test_that("seqEqual", {
    # Ignore gaps
    expect_true(seqEqual("ATG-C", "AT--C"))
    expect_true(seqEqual("ATGGC", "ATGGN"))
    expect_false(seqEqual("AT--T", "ATGGC"))
    
    # Ignore only Ns
    expect_false(seqEqual("ATG-C", "AT--C", ignore="N"))
    expect_true(seqEqual("ATGGC", "ATGGN", ignore="N"))
    expect_false(seqEqual("AT--T", "ATGGC", ignore="N"))
})

test_that("translateDNA", {
    expect_equal(
        translateDNA(df$JUNCTION[1:3]),
        c("CARDRSTPWRRGIASTTVRTSW", "CARDLLWSVLLTGYYSYGMDAW", "CARDLLWSVLLTGYYSYGMDAW"))
    
    expect_equal(
        translateDNA(df$JUNCTION[1:3], trim=TRUE),
        c("ARDRSTPWRRGIASTTVRTS", "ARDLLWSVLLTGYYSYGMDA", "ARDLLWSVLLTGYYSYGMDA"))
    expect_equal(translateDNA("ACTGACTCGA"), "TDS")
})

test_that("maskSeqGaps", {
    expect_equal(maskSeqGaps(c("ATG-C", "CC..C")),
                 c("ATGNC", "CCNNC"))
    
    expect_equal(maskSeqGaps("--ATG-C-"), "NNATGNCN")
    expect_equal(maskSeqGaps("--ATG-C-", outer_only=TRUE), "NNATG-CN")
})


test_that("maskSeqEnds", {
    # Default behavior uniformly masks ragged ends
    seq <- c("CCCCTGGG", "NAACTGGN", "NNNCTGNN")
    expect_equal(maskSeqEnds(seq),c("NNNCTGNN", "NNNCTGNN", "NNNCTGNN"))
    
    # Does nothing
    expect_equal(maskSeqEnds(seq, max_mask=0), c("CCCCTGGG", "NAACTGGN", "NNNCTGNN"))
    
    # Cut ragged sequence ends
    expect_equal(maskSeqEnds(seq, trim=TRUE), c("CTG", "CTG", "CTG"))
    
    # Set max_mask to limit extent of masking and trimming
    maskSeqEnds(seq, max_mask=1)
    expect_equal(maskSeqEnds(seq, max_mask=1, trim=TRUE), c("CCCTGG", "AACTGG", "NNCTGN"))
})

test_that("collapseDuplicates", {
    # Example Change-O data.frame
    df <- data.frame(SEQUENCE_ID=LETTERS[1:4],
                     SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
                     TYPE=c("IgM", "IgG", "IgG", "IgA"),
                     SAMPLE=c("S1", "S1", "S2", "S2"),
                     COUNT=1:4,
                     stringsAsFactors=FALSE)
    
    # Annotations are not parsed if neither text_fields nor num_fields is specified
    # The retained sequence annotations will be random
    obs <- collapseDuplicates(df, verbose=F)
    exp <- data.frame(
        "SEQUENCE_ID" = c("C", "A"),
        "SEQUENCE_IMGT" = c("NAACTGGN", "CCCCTGGG"),
        "TYPE" = c("IgG","IgM"),
        "SAMPLE" = c("S2", "S1"),
        "COUNT" = c(3,1),
        stringsAsFactors = F
    )
    expect_equal(obs, exp)
    
    # Unique text_fields annotations are combined into a single string with ","
    # num_fields annotations are summed
    # Ambiguous duplicates are discarded
    obs <- collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       verbose=F)
    exp$TYPE <- c("IgG","IgG,IgM")
    exp$COUNT <- c(3,3)
    expect_equal(obs, exp)
    
    
    # Use alternate delimiter for collapsing textual annotations
    obs <- collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       sep="/", verbose=F)
    exp$TYPE <- c("IgG","IgG/IgM")
    expect_equal(obs, exp)
    
    # Add count of duplicates
    obs <- collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       add_count=TRUE, verbose=F)
    exp$TYPE <- c("IgG","IgG,IgM")
    exp$COLLAPSE_COUNT <- c(1,2)
    expect_equal(obs, exp)
    
    # Masking ragged ends may impact duplicate removal
    df$SEQUENCE_IMGT <- maskSeqEnds(df$SEQUENCE_IMGT)
    obs <- collapseDuplicates(df, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       add_count=TRUE, verbose=F)    
    exp <- data.frame(
        "SEQUENCE_ID" = "A",
        "SEQUENCE_IMGT" = "NNNCTGNN",
        "TYPE" = "IgA,IgG,IgM",
        "SAMPLE" = "S1,S2",
        "COUNT" = 10,
        "COLLAPSE_COUNT" = 4,
        stringsAsFactors = F
    ) 
    expect_equal(obs, exp)
    
})

test_that("extractVRegion", {
    clone <- subset(df, CLONE == 164)
    
    # Get all regions
    obs <- extractVRegion(clone$SEQUENCE_IMGT)
    expect_equal(dim(unique(obs)), c(7, 5))
    expect_equal(colnames(obs),
                 c("FWR1", "CDR1", "FWR2", "CDR2", "FWR3"))
    
    fwr1 <- c(
        "GAGGTGCAGCTGGTGGAGTCTGG.GGA...GGCTTGGTACAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCT",
        "GAGGTGCAGCTGGTGGAGTCTGG.GGA...GGCTTGGTACAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCT",
        "GAGGTGCAGCTGGTGGAGTCTGG.GGA...GGCTTGGTACAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCT"                     
    )
    expect_equal(obs[1:3,"FWR1"], fwr1)
    
    cdr1 <- c(
        "GGATTCACCTTC............AGTAGTTATGAA",
        "GGATTCACCTTC............AGTAGTTATGAA",
        "GGATTCACCTTC............AGTAGTTATGAA",
        "GGATTCACCTTC............AGTAGTTATGAA"
    )
    expect_equal(obs[4:7,"CDR1"], cdr1)
    
    fwr2 <- c(
        "ATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATAC",
        "ATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATAC",
        "ATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATAC"
    )
    expect_equal(obs[8:10,"FWR2"], fwr2)
    
    cdr2 <- c(
        "ATTAGTAGTAGT......GGTAGTACCATA",
        "ATTAGTAGTAGT......GGTAGTACCATA",
        "ATTAGTAGTAGT......GGTAGTACCATA",
        "ATTAGTAGTAGT......GGTAGTACCATA"
    )
    expect_equal(obs[11:14,"CDR2"], cdr2)
    
    fwr3 <- c(
        "TACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGT",
        "TACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGT",
        "TACTACGCAGACTCTGTGAAG...GGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTTTATTACTGT"
    )
    expect_equal(obs[15:17,"FWR3"], fwr3)
    
    # Get single region
    obs <- extractVRegion(clone$SEQUENCE_IMGT[1:3], "FWR1")
    expect_equal(obs[1:3], fwr1)
    
    # Get all CDRs
    obs <- extractVRegion(clone$SEQUENCE_IMGT, c("CDR1", "CDR2"))
    expect_equal(obs[1:4,"CDR1"],cdr1)
    expect_equal(obs[11:14,"CDR2"],cdr2)
    
    # Get all FWRs
    obs <- extractVRegion(clone$SEQUENCE_IMGT, c("FWR1", "FWR2", "FWR3"))
    expect_equal(obs[1:3,"FWR1"],fwr1)
    expect_equal(obs[8:10,"FWR2"],fwr2)
    expect_equal(obs[15:17,"FWR3"],fwr3)
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