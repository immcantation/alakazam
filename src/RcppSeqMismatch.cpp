#include <Rcpp.h>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <vector>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


inline std::vector<char> parseIgnore(CharacterVector ignore) {
    std::vector<char> ignore_chars;
    ignore_chars.reserve(ignore.length());

    for (int i = 0; i < ignore.length(); i++) {
        if (ignore[i] == NA_STRING) {
            continue;
        }

        std::string ig = as<std::string>(ignore[i]);
        if (ig.length() > 0) {
            ignore_chars.push_back((char)std::toupper((unsigned char)ig[0]));
        }
    }

    return ignore_chars;
}


inline bool isIgnored(char c, const std::vector<char>& ignore) {
    char up = (char)std::toupper((unsigned char)c);

    for (int i = 0; i < (int)ignore.size(); i++) {
        if (up == ignore[i]) {
            return TRUE;
        }
    }

    return FALSE;
}


inline int countMismatches(std::string sample,
                          std::string germline,
                          const std::vector<char>& ignore) {
    int sample_len = sample.length();
    int germline_len = germline.length();
    int len = std::min(sample_len, germline_len);
    int count = 0;

    for (int i = 0; i < len; i++) {
        char sample_char = (char)std::toupper((unsigned char)sample[i]);
        char germline_char = (char)std::toupper((unsigned char)germline[i]);

        if (isIgnored(sample_char, ignore) || isIgnored(germline_char, ignore)) {
            continue;
        }

        if (sample_char != germline_char) {
            count++;
        }
    }

    return count;
}


inline std::vector<int> findMismatchPositions(std::string sample,
                                              std::string germline,
                                              const std::vector<char>& ignore) {
    int sample_len = sample.length();
    int germline_len = germline.length();
    int len = std::min(sample_len, germline_len);

    std::vector<int> pos;

    for (int i = 0; i < len; i++) {
        char sample_char = (char)std::toupper((unsigned char)sample[i]);
        char germline_char = (char)std::toupper((unsigned char)germline[i]);

        if (isIgnored(sample_char, ignore) || isIgnored(germline_char, ignore)) {
            continue;
        }

        if (sample_char != germline_char) {
            pos.push_back(i + 1);
        }
    }

    return pos;
}


//' Count mismatches between sample and germline sequences.
//'
//' \code{seqMismatchCountRcpp} counts Hamming-style mismatches between paired sample and
//' germline sequences, excluding ignored characters.
//'
//' @param    samples    character vector containing sample sequences.
//' @param    germlines  character vector containing germline sequences. If length
//'                      one, the germline is recycled across all samples.
//' @param    ignore     vector of characters to ignore when counting mismatches.
//'                      Default is to ignore c("N", ".", "-").
//'
//' @return   Integer vector of mismatch counts.
//'
//' @details  Comparisons are case-insensitive. Sequences of unequal length are
//'           compared through the length of the shorter sequence.
//'
//' @export
// [[Rcpp::export]]
IntegerVector seqMismatchCountRcpp(CharacterVector samples,
                                CharacterVector germlines,
                                CharacterVector ignore=CharacterVector::create("N", ".", "-")) {
    int n = samples.length();
    int m = germlines.length();

    if (m != 1 && n != m) {
        stop("Number of input sequences does not match number of germlines.");
    }

    std::vector<char> ignore_chars = parseIgnore(ignore);
    IntegerVector counts(n);

    for (int i = 0; i < n; i++) {
        if (samples[i] == NA_STRING) {
            counts[i] = NA_INTEGER;
            continue;
        }

        int germline_i = (m == 1) ? 0 : i;

        if (germlines[germline_i] == NA_STRING) {
            counts[i] = NA_INTEGER;
            continue;
        }

        std::string sample = as<std::string>(samples[i]);
        std::string germline = as<std::string>(germlines[germline_i]);

        counts[i] = countMismatches(sample, germline, ignore_chars);
    }

    return counts;
}


//' Count mismatches between samples and germlines.
//'
//' \code{seqMismatchMatrixRcpp} counts Hamming-style mismatches between each sample
//' and each germline sequence, excluding ignored characters.
//'
//' @param    samples    character vector containing sample sequences.
//' @param    germlines  character vector containing germline sequences.
//' @param    ignore     vector of characters to ignore when counting mismatches.
//'                      Default is to ignore c("N", ".", "-").
//'
//' @return   Integer matrix of mismatch counts, with rows corresponding to
//'           samples and columns corresponding to germlines.
//'
//' @details  Comparisons are case-insensitive. Sequences of unequal length are
//'           compared through the length of the shorter sequence.
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix seqMismatchMatrixRcpp(CharacterVector samples,
                                      CharacterVector germlines,
                                      CharacterVector ignore=CharacterVector::create("N", ".", "-")) {
    int n = samples.length();
    int m = germlines.length();

    std::vector<char> ignore_chars = parseIgnore(ignore);
    IntegerMatrix counts(n, m);

    for (int i = 0; i < n; i++) {
        if (samples[i] == NA_STRING) {
            for (int j = 0; j < m; j++) {
                counts(i, j) = NA_INTEGER;
            }
            continue;
        }

        std::string sample = as<std::string>(samples[i]);

        for (int j = 0; j < m; j++) {
            if (germlines[j] == NA_STRING) {
                counts(i, j) = NA_INTEGER;
                continue;
            }

            std::string germline = as<std::string>(germlines[j]);
            counts(i, j) = countMismatches(sample, germline, ignore_chars);
        }
    }

    Rcpp::List dimnames = Rcpp::List::create(samples.attr("names"),
                                             germlines.attr("names"));
    counts.attr("dimnames") = dimnames;

    return counts;
}


//' Locate mismatches between sample and germline sequences.
//'
//' \code{seqMismatchPositionsRcpp} identifies Hamming-style mismatch positions between
//' paired sample and germline sequences, excluding ignored characters.
//'
//' @param    samples    character vector containing sample sequences.
//' @param    germlines  character vector containing germline sequences. If length
//'                      one, the germline is recycled across all samples.
//' @param    ignore     vector of characters to ignore when locating mismatches.
//'                      Default is to ignore c("N", ".", "-").
//'
//' @return   List of integer vectors containing 1-based mismatch positions.
//'
//' @details  Comparisons are case-insensitive. Sequences of unequal length are
//'           compared through the length of the shorter sequence.
//'
//' @export
// [[Rcpp::export]]
List seqMismatchPositionsRcpp(CharacterVector samples,
                           CharacterVector germlines,
                           CharacterVector ignore=CharacterVector::create("N", ".", "-")) {
    int n = samples.length();
    int m = germlines.length();

    if (m != 1 && n != m) {
        stop("Number of input sequences does not match number of germlines.");
    }

    std::vector<char> ignore_chars = parseIgnore(ignore);
    List positions(n);

    for (int i = 0; i < n; i++) {
        if (samples[i] == NA_STRING) {
            positions[i] = R_NilValue;
            continue;
        }

        int germline_i = (m == 1) ? 0 : i;

        if (germlines[germline_i] == NA_STRING) {
            positions[i] = R_NilValue;
            continue;
        }

        std::string sample = as<std::string>(samples[i]);
        std::string germline = as<std::string>(germlines[germline_i]);

        positions[i] = wrap(findMismatchPositions(sample, germline, ignore_chars));
    }

    return positions;
}
