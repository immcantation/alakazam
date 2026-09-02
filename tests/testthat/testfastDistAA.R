#### fastDistAA vs pairwiseDist ####

test_that("fastDistAA matches pairwiseDist for amino acid sequences", {
  aa_mat <- alakazam::getAAMatrix(gap=0)
  
  # --- 10 random sequences of length 20 contain only 20 standard AA  ---
  AAS  <- c("A","C","D","E","F","G","H","I","K","L",
            "M","N","P","Q","R","S","T","V","W","Y")
  
  seqs <- replicate(10, paste(sample(AAS, 20, replace=TRUE), collapse=""))

  fast_distAA <- alakazam:::fastDistAA(seqs)
  fast_distAA <- as.matrix(fast_distAA) 
  pw_dist<- alakazam::pairwiseDist(seq = seqs, dist_mat = aa_mat)
  
  # Same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
  expect_equal(fast_distAA, pw_dist, check.attributes=F)
  
  
  # --- 6 sequences of length 10 contain only 20 standard AA and X, ., -, *---
  # "X", "-" and "." match everything.  stop codon "*" match itself.
  
  AAS  <- c("A","C","D","E","F","G","H","I","K","L",
            "M","N","P","Q","R","S","T","V","W","Y","X",".","-","*")
  
  seqs <- c("ACDEF*GHIK",
            "AC-XF*GHIK",
            "A.X*FTGHIK",
            "MNPQRSTVWY",
            "*-.XRSTVWY",
            "*CDERSTVWY")
  
  fast_distAA <- alakazam:::fastDistAA(seqs)
  fast_distAA <- as.matrix(fast_distAA) 
  pw_dist<- alakazam::pairwiseDist(seq = seqs, dist_mat = aa_mat)
  
  # Same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
  expect_equal(fast_distAA, pw_dist, check.attributes=F)
  
  
  # --- 1 sequence of length 10 contain only 20 standard AA and X, ., -, *---

  seqs <- c("ACDEF*GHIK")
  fast_distAA <- alakazam:::fastDistAA(seqs)
  fast_distAA <- as.matrix(fast_distAA) 
  pw_dist<- alakazam::pairwiseDist(seq = seqs, dist_mat = aa_mat)
  
  # Same results when comparing to pairwiseDist with check.attributes=F (ignoring dimnames)
  expect_equal(fast_distAA, pw_dist, check.attributes=F)
  
  
  # --- empty sequence list ---
  seq_empty <- c()
  
  expect_error(
    fast_distAA <- alakazam:::fastDistAA(seq_empty),
    "Amino acid sequence list is empty"
  )
  
})