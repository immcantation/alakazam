#### equalityGraph ####

# Random sequences with many ignored characters, so that equal pairs and ambiguous
# sequences are common. Each sequence is one of a few variants of a base sequence,
# with random positions masked.
.randomSeqs <- function(n, L, p_mask, mask=c("N", "-", ".", "?"), lengths=L) {
    base <- sample(c("A", "C", "G", "T"), max(lengths), replace=TRUE)
    variants <- replicate(sample(1:4, 1), {
        v <- base
        m <- runif(length(v)) < 0.3
        v[m] <- sample(c("A", "C", "G", "T"), sum(m), replace=TRUE)
        v
    }, simplify=FALSE)
    vapply(seq_len(n), function(i) {
        len <- lengths[sample.int(length(lengths), 1)]
        s <- variants[[sample.int(length(variants), 1)]][seq_len(len)]
        m <- runif(length(s)) < p_mask
        s[m] <- sample(c(mask, mask[1], mask[1], "n"), sum(m), replace=TRUE)
        paste(s, collapse="")
    }, character(1))
}

test_that("equalityGraph agrees with pairwiseEqual", {
    set.seed(1)
    for (r in 1:50) {
        ignore <- if (r %% 5 == 0) "N" else c("N", "-", ".", "?")
        s <- unique(.randomSeqs(sample(2:60, 1), sample(3:20, 1), runif(1, 0, 0.6),
                                lengths=if (r %% 7 == 0) 8:10 else 10))
        d <- pairwiseEqual(s, ignore=ignore)
        g <- equalityGraph(s, ignore=ignore)
        k <- length(s)

        # Equal pairs, as adjacency lists
        expect_equal(g$n_pairs, (sum(d) - k) / 2)
        expect_identical(g$degree, as.integer(rowSums(d) - 1))
        adj <- lapply(seq_len(k), function(i) {
            if (g$offset[i + 1] > g$offset[i]) { g$adj[(g$offset[i] + 1):g$offset[i + 1]] } else { integer(0) }
        })
        expect_identical(adj, lapply(seq_len(k), function(i) setdiff(which(d[i, ]), i)))

        # Ambiguous: equal to two sequences that are not equal to each other
        clique <- vapply(seq_len(k), function(i) { x <- which(d[i, ]); all(d[x, x]) }, logical(1))
        expect_identical(g$ambiguous, !clique)

        # Components numbered as igraph numbers them
        ig <- igraph::components(igraph::graph_from_adjacency_matrix(d))$membership
        expect_identical(g$component, as.integer(ig))
    }
})

test_that("equalityGraph example", {
    # ANGT is equal to ACGT and AGGT, which are not equal to each other
    g <- equalityGraph(c("ACGT", "ANGT", "AGGT", "ACGTA", "TTTT", "TTTN"))
    expect_equal(g$n_pairs, 3)
    expect_identical(g$degree, c(1L, 2L, 1L, 0L, 1L, 1L))
    expect_identical(g$ambiguous, c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
    expect_identical(g$set, c(1L, NA, 3L, 4L, 5L, 5L))
    expect_identical(g$component, c(1L, 1L, 1L, 2L, 3L, 3L))
})


#### collapseDuplicates against the dense implementation ####

# .collapseDuplicatesDense, in helper-collapseDuplicatesDense.R, is the implementation
# that used the dense pairwiseEqual matrix. Results must be identical, including for
# dry=TRUE, ties, ambiguous sequences, and the case where all or all but one
# sequence is ambiguous.
test_that("collapseDuplicates is identical to the dense implementation", {
    set.seed(2)
    n_amb_all <- 0
    for (r in 1:400) {
        n <- sample(2:30, 1)
        seqs <- .randomSeqs(n, sample(3:10, 1),
                            if (r > 250) runif(1, 0.5, 0.9) else runif(1, 0, 0.5))
        if (r %% 10 == 0) { seqs[sample(n, 1)] <- strrep("N", nchar(seqs[1])) }
        db <- data.frame(sequence_id=sprintf("S%02d", 1:n), sequence_alignment=seqs,
                         c_call=sample(c("IGHM", "IGHG1,IGHA1", NA), n, replace=TRUE),
                         duplicate_count=sample(c(1:5, NA), n, replace=TRUE),
                         consensus_count=as.numeric(sample(1:5, n, replace=TRUE)),
                         germline=seqs, stringsAsFactors=FALSE)
        for (dry in c(FALSE, TRUE)) {
            args <- list(db, id="sequence_id", seq="sequence_alignment",
                         text_fields=if (r %% 3 == 0) "c_call" else NULL,
                         num_fields=c("duplicate_count", "consensus_count"),
                         seq_fields=if (r %% 5 == 0) "germline" else NULL,
                         add_count=r %% 2 == 0, ignore=c("N", "-", ".", "?"), sep=",",
                         dry=dry, verbose=FALSE)
            expect_identical(suppressWarnings(do.call(.collapseDuplicatesSingle, args)),
                             suppressWarnings(do.call(.collapseDuplicatesDense, args)))
        }
        if (sum(!equalityGraph(unique(seqs))$ambiguous[match(seqs, unique(seqs))]) <= 1) {
            n_amb_all <- n_amb_all + 1
        }
    }
    # The all-ambiguous case was exercised
    expect_gt(n_amb_all, 10)
})
