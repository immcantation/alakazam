#include <Rcpp.h>
using namespace Rcpp;

inline uint8_t code_char(char c){
  switch(c){
    case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3;
    case 'N': return 4; case '?': return 5;
    default:  return 255;
  }
}

// [[Rcpp::export]]
int countSeqsWithInvalidBases_rcpp(CharacterVector seqs) {
  int N = seqs.size();
  int bad = 0;

  for (int i = 0; i < N; ++i) {
    if (seqs[i] == NA_STRING) {
      bad++;
      continue;
    }

    std::string s = as<std::string>(seqs[i]);
    bool invalid = false;
    for (char ch : s) {
      char up = (char)std::toupper((unsigned char)ch);
      if (code_char(up) == 255) {
        invalid = true;
        break;
      }
    }

    if (invalid) bad++;
  }

  return bad;
}

// [[Rcpp::export]]
IntegerVector fastDist_rcpp(CharacterVector seqs) {
  int N = seqs.size();

  if (N == 0) stop("empty input");

  std::string s0 = as<std::string>(seqs[0]);
  int L = (int)s0.size();
  for (int i = 0; i < N; ++i) {
    if ((int)std::string(as<std::string>(seqs[i])).size() != L)
      stop("All sequences must have the same length");
  }

  // encode to uint8: N x L, row-major in a flat vector
  std::vector<uint8_t> enc((size_t)N * L);
  for (int i = 0; i < N; ++i) {
    std::string s = as<std::string>(seqs[i]);
    for (int p = 0; p < L; ++p) {
      uint8_t c = code_char(s[p]);
      if (c == 255) stop("Only A,C,G,T,N,? are allowed");
      enc[(size_t)i * L + p] = c;
    }
  }

  // lower-triangle vector (column-major, matching R's dist storage):
  // element at row i, col j (i > j, 0-indexed) maps to position:
  //   j*N - j*(j+1)/2 + (i-j) - 1
  size_t tri_size = (size_t)N * (N - 1) / 2;
  IntegerVector tri(tri_size); // zero-initialized

  auto tri_idx = [&](int i, int j) -> size_t {
    if (i < j) std::swap(i, j);
    return (size_t)j * N - (size_t)j * (j + 1) / 2 + (i - j) - 1;
  };

  for (int p = 0; p < L; ++p) {
    std::vector<int> A, Cb, G, Tt, Ns, Q;
    A.reserve(N); Cb.reserve(N); G.reserve(N); Tt.reserve(N); Ns.reserve(N); Q.reserve(N);

    // bucket row indices by symbol at column p
    for (int i = 0; i < N; ++i) {
      switch (enc[(size_t)i * L + p]) {
        case 0: A.push_back(i);  break;
        case 1: Cb.push_back(i); break;
        case 2: G.push_back(i);  break;
        case 3: Tt.push_back(i); break;
        case 4: Ns.push_back(i); break;
        case 5: Q.push_back(i);  break;
      }
    }

    // increment lower-triangle match counts for pairs within the same bucket
    auto bump_within = [&](const std::vector<int>& v) {
      int m = (int)v.size();
      for (int ii = 1; ii < m; ++ii)
        for (int jj = 0; jj < ii; ++jj)
          tri[tri_idx(v[ii], v[jj])] += 1;
    };

    // increment lower-triangle match counts for all cross-bucket pairs
    auto bump_pairs = [&](const std::vector<int>& X, const std::vector<int>& Y) {
      for (int xi : X)
        for (int yj : Y)
          if (xi != yj)
            tri[tri_idx(xi, yj)] += 1;
    };

    // same known base
    bump_within(A);  bump_within(Cb);  bump_within(G);  bump_within(Tt);

    // N with known (both directions), and N with N (N matches any base)
    std::vector<int> K; K.reserve(A.size()+Cb.size()+G.size()+Tt.size());
    K.insert(K.end(), A.begin(),  A.end());
    K.insert(K.end(), Cb.begin(), Cb.end());
    K.insert(K.end(), G.begin(),  G.end());
    K.insert(K.end(), Tt.begin(), Tt.end());
    bump_pairs(Ns, K);
    // bump_pairs(K, Ns); //  The bump_pairs(Ns, K) call already handles all N-vs-known-base pairs in both directions because tri_idx normalizes (i,j) and (j,i) to the same lower-triangle position.
    bump_within(Ns);  // N-N is a match (consistent with getDNAMatrix(gap=0))

    // ? with ?
    bump_within(Q);
  }

  // convert matches -> distances
  for (int k = 0; k < tri_size; ++k)
    tri[k] = L - tri[k];

  return tri;
}
