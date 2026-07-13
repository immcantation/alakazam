#### fastDist vs pairwiseDist ####

test_that("fastDist matches pairwiseDist for ATCG sequences", {
    dna_mat <- alakazam::getDNAMatrix(gap=0)

    # --- Basic ATCG cases ---
    seqs <- c(
        "ACGTACGT",  # seq1
        "ACGTACGT",  # seq2: identical to seq1 (0 mismatches)
        "ACGTACGC",  # seq3: 1 mismatch from seq1 (last position T->C)
        "TTTTTTTT"   # seq4: 6 mismatches from seq1 (positions 1,2,3,5,6,7)
    )

    fast_dist <- alakazam:::fastDist(seqs)
    fast_dist <- as.matrix(fast_dist) 
    pw_dist     <- alakazam::pairwiseDist(seqs, dna_mat)

    # Same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
    expect_equal(fast_dist, pw_dist, check.attributes=F)

    # --- N behaviour ---
    # N represents any nucleotide. Distance 0 vs any known base

    seqs_n <- c("ACGN", "ACGN", "ACGA", "ACGT")
    fast_n <- alakazam::fastDist(seqs_n)
    fast_n <- as.matrix(fast_n)
    pw_n   <- alakazam::pairwiseDist(seqs_n, dna_mat)

    # Same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
    expect_equal(fast_n, pw_n, check.attributes=F)

    # --- ? behaviour ---
    # ? means missing data: matches only itself, mismatches everything else

    seqs_q <- c("ACG?", "ACG?", "ACGA", "ACGN")
    fast_q <- alakazam::fastDist(seqs_q)
    fast_q <- as.matrix(fast_q)
    pw_q   <- alakazam::pairwiseDist(seqs_q, dna_mat)

    # Same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
    expect_equal(fast_q, pw_q, check.attributes=F)

    # --- Mixed N, ?, and ATCG: full matrix matches pairwiseDist ---
    seqs_mixed <- c("ACGTNACGT?", "ACGTNACGT?", "ACGTAACGTA", "TTTTTTTTTT")
    fast_mixed <- alakazam::fastDist(seqs_mixed)
    fast_mixed <- as.matrix(fast_mixed)
    pw_mixed   <- alakazam::pairwiseDist(seqs_mixed, dna_mat)

    # Expect same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
    expect_equal(fast_mixed, pw_mixed, check.attributes=FALSE)

    # --- Single sequence: 1x1 matrix, diagonal = 0 ---
    fast_single <- alakazam::fastDist("ACGT")
    fast_single <- as.matrix(fast_single)
    single <- alakazam::pairwiseDist("ACGT", dna_mat)

    # Expect same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
    expect_equal(fast_single, single, check.attributes = FALSE)
})