ExampleDb <- file.path("..", "data-tests", "ExampleDb.gz")
db <- readChangeoDb(ExampleDb)

# Load test database
e1 <- new.env()
#load(file.path("tests", "data-tests", "TestDb.rda"), envir=e1)
load(file.path("..", "data-tests", "TestDb.rda"), envir=e1)
db2 <- get("TestDb", envir=e1)
rm(e1)

test_that("seqDist: short toy sequences", {
    
    expect_equal(seqDist("AC-A", "AC-G", getDNAMatrix(gap=0)), 1)
    expect_equal(seqDist("AC--A", "ACATG", getDNAMatrix(gap=0)), 1)
    expect_equal(seqDist("AC--A", "ACATG", getDNAMatrix(gap=1)), 3)
    expect_equal(seqDist("AC--A", "ACATG", getDNAMatrix(gap=-1)), 2)
    
    # AC--AAC--A
    # ACATGACATG
    # **--.**--.
    expect_equal(seqDist("AC--AAC--A", "ACATGACATG", getDNAMatrix(gap=0)), 2)
    expect_equal(seqDist("AC--AAC--A", "ACATGACATG", getDNAMatrix(gap=1)), 6)
    expect_equal(seqDist("AC--AAC--A", "ACATGACATG", getDNAMatrix(gap=-1)), 4)
    
    
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

test_that("seqDist: long IMGT-gapped sequences", {
    # Define test sequence set
    seq <- substr(db2$SEQUENCE_IMGT[1:4], 1, 312)
    germ <- substr(db2$GERMLINE_IMGT_D_MASK[1:4], 1, 312)
    # Replace dots with Ns
    seq_n <- gsub(".", "N", seq, fixed=TRUE)
    germ_n <- gsub(".", "N", germ, fixed=TRUE)
    
    # Region ranges
    regions <- list("SEQ"=c(1, 312),
                    "FWR1"=c(79, 114),
                    "CDR1"=c(115, 165),
                    "FWR2"=c(115, 165),
                    "CDR2"=c(166, 195),
                    "FWR3"=c(196, 312))
    
    # Expected mutations
    expected <- list("SEQ"=c(13, 19, 19, 55),
                     "FWR1"=c(2, 6, 6, 4),
                     "CDR1"=c(1, 5, 5, 3),
                     "FWR2"=c(2, 1, 1, 3),
                     "CDR2"=c(3, 3, 3, 6),
                     "FWR3"=c(5, 4, 4, 39))
    
    # Test full V region sequence
    #cat("Full V-region (", paste(regions[["SEQ"]], collapse=":"), ") distance:\n", sep="")
    
    #cat("  With dots:\n")
    d <- mapply(seqDist, seq, germ, MoreArgs=list(dist_mat=getDNAMatrix(gap=0)), 
                USE.NAMES=FALSE)
    #d <- mapply(function(x, y) { sum(seqinr::s2c(x) != seqinr::s2c(y)) }, seq, germ, USE.NAMES=FALSE)
    expect_equal(d, expected[["SEQ"]], info="Full V-region with dots")
    
    #cat("  With Ns:\n")
    d <- mapply(seqDist, seq_n, germ_n, MoreArgs=list(dist_mat=getDNAMatrix(gap=0)), 
                USE.NAMES=FALSE)
    #d <- mapply(function(x, y) { sum(seqinr::s2c(x) != seqinr::s2c(y)) }, seq, germ, USE.NAMES=FALSE)
    expect_equal(d, expected[["SEQ"]], info="Full V-region with Ns")
    
    # Test by region
    for (n in c("FWR1", "CDR1", "FWR2", "CDR2", "FWR3")) {
        #cat(n, " (", paste(regions[[n]], collapse=":"), ") distance:\n", sep="")
        
        # Define substrings
        seq_sub <- extractVRegion(seq, n)
        germ_sub <- extractVRegion(germ, n)
        seq_n_sub <- extractVRegion(seq_n, n)
        germ_n_sub <- extractVRegion(germ_n, n)
        
        #cat("  With dots:\n")
        d <- mapply(seqDist, seq_sub, germ_sub, MoreArgs=list(dist_mat=getDNAMatrix(gap=0)), 
                    USE.NAMES=FALSE)
        #d <- mapply(function(x, y) { sum(seqinr::s2c(x) != seqinr::s2c(y)) }, seq_sub, germ_sub, USE.NAMES=FALSE)
        expect_equal(d, expected[[n]], info=paste(n, "with dots"))
        
        #cat("  With Ns:\n")
        d <- mapply(seqDist, seq_n_sub, germ_n_sub, MoreArgs=list(dist_mat=getDNAMatrix(gap=0)), 
                    USE.NAMES=FALSE)
        #d <- mapply(function(x, y) { sum(seqinr::s2c(x) != seqinr::s2c(y)) }, seq_sub, germ_sub, USE.NAMES=FALSE)
        expect_equal(d, expected[[n]], info=paste(n, "with Ns"))
    }
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
        translateDNA(db$JUNCTION[1:3]),
        c("CARDRSTPWRRGIASTTVRTSW", "CARDLLWSVLLTGYYSYGMDAW", "CARDLLWSVLLTGYYSYGMDAW"))
    
    expect_equal(
        translateDNA(db$JUNCTION[1:3], trim=TRUE),
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
    db <- data.frame(SEQUENCE_ID=LETTERS[1:5],
                     SEQUENCE_IMGT=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN", "NAACTGNG"),
                     TYPE=c("IgM", "IgG", "IgG", "IgA", "IgG"),
                     SAMPLE=c("S1", "S1", "S2", "S2","S2"),
                     COUNT=1:5,
                     stringsAsFactors=FALSE)
    
    # Annotations are not parsed if neither text_fields nor num_fields is specified
    # The retained sequence annotations will be random
    obs <- collapseDuplicates(db, verbose=F)
    exp <- data.frame(
        "SEQUENCE_ID" = c("C", "A"),
        "SEQUENCE_IMGT" = c("NAACTGGN", "CCCCTGGG"),
        "TYPE" = c("IgG","IgM"),
        "SAMPLE" = c("S2", "S1"),
        "COUNT" = c(3,1),
        stringsAsFactors = F
    )
    expect_equivalent(obs, exp[2:1,])
    
    obs_dry <- collapseDuplicates(db[-5,], verbose=F, dry=T)
    expect_equal(obs_dry$COLLAPSE_CLASS, c("duplicated", "duplicated", "unique", "ambiguous"))
    
    expect_equal(sort(obs_dry[obs_dry$COLLAPSE_PASS,"SEQUENCE_ID"]), 
                 sort(obs$SEQUENCE_ID))
    
    ## Try messing up order
    ## C comes first
    obs_dry <- collapseDuplicates(db, verbose=F, dry=T)
    expect_equal(sort(obs_dry[obs_dry$COLLAPSE_PASS,"SEQUENCE_ID"]), 
                 sort(obs$SEQUENCE_ID))
    ## E comes first
    obs_dry <- collapseDuplicates(db[nrow(db):1,], verbose=F, dry=T)
    expect_equal(sort(obs_dry[obs_dry$COLLAPSE_PASS,"SEQUENCE_ID"]), 
                 c("A","E"))
    
    # Unique text_fields annotations are combined into a single string with ","
    # num_fields annotations are summed
    # Ambiguous duplicates are discarded
    obs <- collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       verbose=F)
    exp$TYPE <- c("IgG","IgG,IgM")
    exp$COUNT <- c(3,3)
    expect_equal(obs, exp)
    
    
    # Use alternate delimiter for collapsing textual annotations
    obs <- collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       sep="/", verbose=F)
    exp$TYPE <- c("IgG","IgG/IgM")
    expect_equal(obs, exp)
    
    # Add count of duplicates
    obs <- collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
                       add_count=TRUE, verbose=F)
    exp$TYPE <- c("IgG","IgG,IgM")
    exp$COLLAPSE_COUNT <- c(1,2)
    expect_equal(obs, exp)
    
    # Masking ragged ends may impact duplicate removal
    db$SEQUENCE_IMGT <- maskSeqEnds(db$SEQUENCE_IMGT)
    obs <- collapseDuplicates(db, text_fields=c("TYPE", "SAMPLE"), num_fields="COUNT", 
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
    clone <- subset(db, CLONE == 164)
    
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