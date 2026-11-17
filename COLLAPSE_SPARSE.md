# collapseDuplicates without the k × k matrix

Branch `collapse_sparse`, based on `fix_oom` @ `5cd70e6`.

## The problem

`.collapseDuplicatesSingle` built `pairwiseEqual(unique(seq))`, a dense logical
matrix of 4 bytes per pair of unique sequences. On deep bulk BCR repertoires a single
group can hold hundreds of thousands of unique sequences: 186 GB of matrix at
k = 223,332 and 973 GB at k = 511,028. Even with the `fix_oom` change, sample 4A6
needed 2.05 TB and was killed at the 23 h limit in `COLLAPSE_DUPLICATES`.

The matrix is almost entirely `FALSE`. Unique sequences can only be equal through
ignored characters, and measured densities are 5 × 10⁻⁵ to 6 × 10⁻⁴.

## The change

- **`src/RcppCollapse.cpp`** (new): internal `equalityGraph()`. It stores only the equal
  pairs, as sorted lists of each sequence's matches, and computes from them everything
  `.collapseDuplicatesSingle` used the matrix for:
  - which sequences are ambiguous (a sequence's matches are not all equal to each other)
  - the duplicate sets
  - connected components, for the all-ambiguous case

  Equality is exactly `seqEqual`'s, compared 8 bytes at a time with early exit.
  Memory grows with the number of equal pairs instead of k².
- **`R/Sequence.R`**: `.collapseDuplicatesSingle` rewritten on `equalityGraph()`,
  covering both `dry = FALSE` and `dry = TRUE`. Per-set `%in%` scans are replaced by
  vectorised indexing. `pairwiseEqual` is unchanged and still exported.
- **Tests**: `tests/testthat/testCollapseDuplicates.R` checks `equalityGraph` against
  `pairwiseEqual` and igraph. It also checks 400 random groups, with `dry` on and off,
  against the old implementation, kept as `helper-collapseDuplicatesDense.R`.
- **`NEWS.md`** entry.

## Results are unchanged

- 18,000 / 18,000 random-group comparisons give `identical()` results against `5cd70e6`.
  These cover `dry` on and off, `fields`, `text_fields` / `seq_fields`, ties, `NA`
  counts, and 617 groups on the all-ambiguous path.
- On six real groups from sample 3E5 (k = 5,011 to 141,758, 528,903 rows),
  `enchantr::findDuplicates` output is identical in all 69 columns.
- The existing alakazam test suite passes.

## Memory and run time

`findDuplicates`, `nproc = 1`, one collapse group of sample 3E5 per row; both versions
run in the same job:

| unique seqs (k) | rows | dense: time | dense: peak heap | sparse: time | sparse: peak heap |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 5,011 | 5,465 | 13 s | 2.5 GB | 3 s | 0.7 GB |
| 19,429 | 20,711 | 212 s | 5.9 GB | 4 s | 0.7 GB |
| 39,285 | 47,542 | 700 s | 12.4 GB | 18 s | 0.9 GB |
| 58,049 | 77,911 | 1,328 s | 26.0 GB | 36 s | 1.0 GB |
| 100,385 | 134,100 | 4,318 s | 76.1 GB | 150 s | 1.2 GB |
| 141,758 | 243,370 | 10,317 s | 151 GB | 298 s | 1.8 GB |

The whole airrflow `COLLAPSE_DUPLICATES` task for sample 4A6 (4,284,341 rows, largest
group k = 511,028), in the `airrflow-5.2.0test` image:

| | dense (`fix_oom`) | sparse |
| --- | --- | --- |
| outcome | killed at the 23 h limit | completed |
| elapsed | > 23 h | 3 h 38 m |
| peak memory (MaxRSS) | 2.05 TB | 32 GB |

Pair comparisons are still O(k²) in time, but cheap: from k = 100,385 to 141,758 the
sparse version's run time grew by the k² factor (1.99×), and the dense version's by
2.39×.
