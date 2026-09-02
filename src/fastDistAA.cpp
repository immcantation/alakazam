#include <Rcpp.h>
using namespace Rcpp;

// A=0 C=1 D=2 E=3 F=4 G=5 H=6 I=7 K=8 L=9
// M=10 N=11 P=12 Q=13 R=14 S=15 T=16 V=17 W=18 Y=19
// X=20 .=21 -=22 (universal wildcards) *=23 (matches only X/./-/itself)
inline uint8_t code_aa(char c){
  switch(c){
    case 'A': return 0;  case 'C': return 1;  case 'D': return 2;
    case 'E': return 3;  case 'F': return 4;  case 'G': return 5;
    case 'H': return 6;  case 'I': return 7;  case 'K': return 8;
    case 'L': return 9;  case 'M': return 10; case 'N': return 11;
    case 'P': return 12; case 'Q': return 13; case 'R': return 14;
    case 'S': return 15; case 'T': return 16; case 'V': return 17;
    case 'W': return 18; case 'Y': return 19;
    case 'X': return 20; case '.': return 21;
    case '-': return 22; case '*': return 23;
    default:  return 255;
  }
}

// [[Rcpp::export]]
int countAASeqsWithInvalidChars_rcpp(CharacterVector seqs) {
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
      if (code_aa(up) == 255) {
        invalid = true;
        break;
      }
    }

    if (invalid) bad++;
  }

  return bad;
}

// [[Rcpp::export]]
IntegerVector fastDistAA_rcpp(CharacterVector seqs) {
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
      char up = (char)std::toupper((unsigned char)s[p]);
      uint8_t c = code_aa(up);
      if (c == 255) stop("Only the 20 standard AAs plus X, ., -, and * are allowed");
      enc[(size_t)i * L + p] = c;
    }
  }

  // lower-triangle vector (column-major, matching R's dist storage):
  // element at row i, col j (i > j, 0-indexed) maps to:
  //   j*N - j*(j+1)/2 + (i-j) - 1
  size_t tri_size = (size_t)N * (N - 1) / 2;
  IntegerVector tri(tri_size); // zero-initialized

  auto tri_idx = [&](int i, int j) -> size_t {
    if (i < j) std::swap(i, j);
    return (size_t)j * N - (size_t)j * (j + 1) / 2 + (i - j) - 1;
  };

  for (int p = 0; p < L; ++p) {
    // one bucket per code (0-19 known AAs, 20=X, 21=., 22=-, 23=*)
    std::vector<int> buckets[24];
    for (int b = 0; b < 24; ++b) buckets[b].reserve(8);

    for (int i = 0; i < N; ++i)
      buckets[enc[(size_t)i * L + p]].push_back(i);

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

    // same known AA
    for (int b = 0; b < 20; ++b)
      bump_within(buckets[b]);

    // X, ., - are universal wildcards: match each other and * (but AA-vs-AA
    // is already covered above, and AA-vs-* never matches)
    std::vector<int> W; // X, ., -
    W.reserve(N);
    for (int b = 20; b <= 22; ++b)
      W.insert(W.end(), buckets[b].begin(), buckets[b].end());

    std::vector<int> C = W; // X, ., -, *
    C.insert(C.end(), buckets[23].begin(), buckets[23].end());

    // any pair drawn from {X, ., -, *} matches, regardless of which symbol
    bump_within(C);

    // any known AA matches X, ., - (but not *)
    std::vector<int> K;
    K.reserve(N);
    for (int b = 0; b < 20; ++b)
      K.insert(K.end(), buckets[b].begin(), buckets[b].end());
    bump_pairs(K, W);
  }

  // convert matches -> distances
  for (size_t k = 0; k < tri_size; ++k)
    tri[k] = L - tri[k];

  return tri;
}