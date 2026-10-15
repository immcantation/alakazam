// Equality graph used by collapseDuplicates
//
// .collapseDuplicatesSingle previously built pairwiseEqual(unique(seq)), a dense
// k x k logical matrix of 4k^2 bytes, and used it only to find, for each sequence,
// the other sequences it is equal to. Distinct sequences can only be equal through
// ignored characters, so the matrix is almost entirely FALSE (density ~1e-4 on bulk
// BCR repertoires) and at k ~ 2e5 to 5e5 it needs hundreds of GB to TB.
//
// equalityGraph() stores only the equal pairs, as sorted adjacency lists, and
// computes from them everything .collapseDuplicatesSingle needs. Memory is
// proportional to the number of equal pairs. Equality is identical to seqEqual().

#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

// Sequences of one length, packed so that a pair can be compared 8 bytes at a time.
// bytes holds each sequence zero-padded to a multiple of 8 bytes; known is 0xFF at
// positions that are not ignore characters and 0x00 at those that are. A pair is
// unequal iff some word has (a ^ b) & known_a & known_b != 0. Padding is zero in
// both arrays, so it never contributes.
struct PackedSeqs {
    std::size_t words = 0;
    std::vector<uint64_t> bytes, known;
};

inline bool packedEqual(const PackedSeqs& p, std::size_t a, std::size_t b) {
    const uint64_t* xa = &p.bytes[a * p.words];
    const uint64_t* xb = &p.bytes[b * p.words];
    const uint64_t* ka = &p.known[a * p.words];
    const uint64_t* kb = &p.known[b * p.words];
    for (std::size_t w = 0; w < p.words; ++w) {
        if ((xa[w] ^ xb[w]) & ka[w] & kb[w]) { return false; }
    }
    return true;
}

// Is the sorted range a a subset of the sorted range b?
inline bool isSubset(const int* a, const int* a_end, const int* b, const int* b_end) {
    while (a != a_end) {
        while (b != b_end && *b < *a) { ++b; }
        if (b == b_end || *b != *a) { return false; }
        ++a; ++b;
    }
    return true;
}

int findRoot(std::vector<int>& parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

}  // namespace


// Equality graph of a set of sequences
//
// Internal. Used by .collapseDuplicatesSingle in place of pairwiseEqual.
//
// @param   seq     character vector of sequences, normally unique(data[[seq]]).
// @param   ignore  vector of characters to ignore when testing for equality. As in
//                  seqEqual, only the first character of each value is used.
//
// @return  A list with
//   n_pairs:    number of unordered pairs i < j that are equal.
//   degree:     integer, the number of other sequences each sequence is equal to.
//   offset, adj: the equal pairs as compressed sparse rows: the sequences equal to
//               sequence i (1-based) are adj[(offset[i] + 1):offset[i + 1]], 1-based
//               and sorted.
//   ambiguous:  logical, TRUE if the sequence is equal to two sequences that are
//               not equal to each other.
//   set:        integer, for sequences that are not ambiguous, the index of the
//               first sequence in its set of duplicates once ambiguous sequences
//               are excluded. NA for ambiguous sequences.
//   component:  integer, connected component of the graph including ambiguous
//               sequences, numbered 1, 2, ... in order of each component's first
//               sequence, as igraph::components numbers them.
// [[Rcpp::export]]
List equalityGraph(CharacterVector seq,
                   CharacterVector ignore=CharacterVector::create("N","-",".","?")) {
    const int k = seq.size();

    bool is_ignore[256] = {false};
    for (R_xlen_t j = 0; j < ignore.size(); j++) {
        const char* s = CHAR(STRING_ELT(ignore, j));
        is_ignore[(unsigned char)s[0]] = true;
    }

    // Sequences of different lengths are never equal, so compare within lengths.
    // seqEqual receives NA as the string "NA", so do the same.
    std::vector<const char*> str(k);
    std::vector<std::size_t> len(k);
    for (int i = 0; i < k; i++) {
        SEXP e = STRING_ELT(seq, i);
        str[i] = (e == NA_STRING) ? "NA" : CHAR(e);
        len[i] = std::strlen(str[i]);
    }
    std::vector<int> by_len(k);
    std::iota(by_len.begin(), by_len.end(), 0);
    std::stable_sort(by_len.begin(), by_len.end(),
                     [&](int a, int b) { return len[a] < len[b]; });

    // Find all equal pairs i < j
    std::vector<std::pair<int, int>> edges;
    for (std::size_t lo = 0; lo < (std::size_t)k; ) {
        std::size_t hi = lo;
        while (hi < (std::size_t)k && len[by_len[hi]] == len[by_len[lo]]) { hi++; }
        const std::size_t n = hi - lo;
        const std::size_t L = len[by_len[lo]];

        PackedSeqs p;
        p.words = (L + 7) / 8;
        p.bytes.assign(n * p.words, 0);
        p.known.assign(n * p.words, 0);
        for (std::size_t s = 0; s < n; s++) {
            const char* c = str[by_len[lo + s]];
            unsigned char* b = reinterpret_cast<unsigned char*>(&p.bytes[s * p.words]);
            unsigned char* m = reinterpret_cast<unsigned char*>(&p.known[s * p.words]);
            for (std::size_t q = 0; q < L; q++) {
                b[q] = (unsigned char)c[q];
                m[q] = is_ignore[(unsigned char)c[q]] ? 0x00 : 0xFF;
            }
        }
        for (std::size_t a = 0; a < n; a++) {
            if (a % 1024 == 0) { Rcpp::checkUserInterrupt(); }
            for (std::size_t b = a + 1; b < n; b++) {
                if (packedEqual(p, a, b)) {
                    int ia = by_len[lo + a];
                    int ib = by_len[lo + b];
                    edges.emplace_back(std::min(ia, ib), std::max(ia, ib));
                }
            }
        }
        lo = hi;
    }

    // Adjacency lists in both directions, sorted
    std::vector<int> deg(k, 0);
    for (const auto& e : edges) { deg[e.first]++; deg[e.second]++; }
    std::vector<std::size_t> off(k + 1, 0);
    for (int i = 0; i < k; i++) { off[i + 1] = off[i] + deg[i]; }
    std::vector<int> adj(off[k]);
    {
        std::vector<std::size_t> pos(off.begin(), off.end() - 1);
        for (const auto& e : edges) {
            adj[pos[e.first]++] = e.second;
            adj[pos[e.second]++] = e.first;
        }
    }
    const double n_pairs = (double)edges.size();
    std::vector<std::pair<int, int>>().swap(edges);
    for (int i = 0; i < k; i++) { std::sort(adj.begin() + off[i], adj.begin() + off[i + 1]); }

    // Closed neighbourhood {i} + adj(i), sorted
    auto closed = [&](int i, std::vector<int>& out) {
        out.assign(adj.begin() + off[i], adj.begin() + off[i + 1]);
        out.insert(std::upper_bound(out.begin(), out.end(), i), i);
    };

    // A sequence is ambiguous if its closed neighbourhood is not a clique, which is
    // the case iff it is not a subset of the closed neighbourhood of every neighbour.
    LogicalVector ambiguous(k, false);
    {
        std::vector<int> ci, cj;
        for (int i = 0; i < k; i++) {
            if (deg[i] == 0) { continue; }
            closed(i, ci);
            for (std::size_t q = off[i]; q < off[i + 1]; q++) {
                int j = adj[q];
                if (deg[j] < deg[i]) { ambiguous[i] = true; break; }
                closed(j, cj);
                if (!isSubset(ci.data(), ci.data() + ci.size(),
                              cj.data(), cj.data() + cj.size())) {
                    ambiguous[i] = true;
                    break;
                }
            }
        }
    }

    // Sets of duplicates. Two equal sequences that are not ambiguous have cliques as
    // closed neighbourhoods, each containing the other, so the neighbourhoods are
    // identical. Each set is therefore well defined and labelled by its first member.
    IntegerVector set(k, NA_INTEGER);
    for (int i = 0; i < k; i++) {
        if (ambiguous[i]) { continue; }
        int m = i;
        for (std::size_t q = off[i]; q < off[i + 1]; q++) {
            if (!ambiguous[adj[q]] && adj[q] < m) { m = adj[q]; }
        }
        set[i] = m + 1;
    }

    // Connected components. Roots are always the smallest member, so numbering roots
    // in order of first appearance numbers components by their first member.
    std::vector<int> parent(k);
    std::iota(parent.begin(), parent.end(), 0);
    for (int i = 0; i < k; i++) {
        for (std::size_t q = off[i]; q < off[i + 1]; q++) {
            int a = findRoot(parent, i);
            int b = findRoot(parent, adj[q]);
            if (a != b) { parent[std::max(a, b)] = std::min(a, b); }
        }
    }
    IntegerVector component(k);
    std::vector<int> label(k, 0);
    int n_comp = 0;
    for (int i = 0; i < k; i++) {
        int r = findRoot(parent, i);
        if (!label[r]) { label[r] = ++n_comp; }
        component[i] = label[r];
    }

    // 1-based adjacency for R
    NumericVector offset_r(k + 1);
    for (int i = 0; i <= k; i++) { offset_r[i] = (double)off[i]; }
    IntegerVector adj_r(adj.size());
    for (std::size_t q = 0; q < adj.size(); q++) { adj_r[q] = adj[q] + 1; }

    return List::create(_["n_pairs"]=n_pairs,
                        _["degree"]=IntegerVector(deg.begin(), deg.end()),
                        _["offset"]=offset_r,
                        _["adj"]=adj_r,
                        _["ambiguous"]=ambiguous,
                        _["set"]=set,
                        _["component"]=component);
}
