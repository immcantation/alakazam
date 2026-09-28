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
                          const std::vector<char>& ignore,
                          bool count_trailing) {
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

    // Sample positions past the end of the germline cannot match it
    if (count_trailing) {
        for (int i = len; i < sample_len; i++) {
            if (!isIgnored(sample[i], ignore)) {
                count++;
            }
        }
    }

    return count;
}


inline std::vector<int> findMismatchPositions(std::string sample,
                                              std::string germline,
                                              const std::vector<char>& ignore,
                                              bool count_trailing) {
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

    // Sample positions past the end of the germline cannot match it
    if (count_trailing) {
        for (int i = len; i < sample_len; i++) {
            if (!isIgnored(sample[i], ignore)) {
                pos.push_back(i + 1);
            }
        }
    }

    return pos;
}


// Internal; the documented, exported interface is seqMismatchCount() in R/Sequence.R.
// [[Rcpp::export]]
IntegerVector seqMismatchCountRcpp(CharacterVector samples,
                                CharacterVector germlines,
                                CharacterVector ignore=CharacterVector::create("N", "-", ".", "?"),
                                bool count_trailing=false) {
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

        counts[i] = countMismatches(sample, germline, ignore_chars, count_trailing);
    }

    return counts;
}


// Internal; the documented, exported interface is seqMismatchMatrix() in R/Sequence.R.
// [[Rcpp::export]]
IntegerMatrix seqMismatchMatrixRcpp(CharacterVector samples,
                                      CharacterVector germlines,
                                      CharacterVector ignore=CharacterVector::create("N", "-", ".", "?"),
                                bool count_trailing=false) {
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
            counts(i, j) = countMismatches(sample, germline, ignore_chars, count_trailing);
        }
    }

    Rcpp::List dimnames = Rcpp::List::create(samples.attr("names"),
                                             germlines.attr("names"));
    counts.attr("dimnames") = dimnames;

    return counts;
}


// Internal; the documented, exported interface is seqMismatchPositions() in R/Sequence.R.
// [[Rcpp::export]]
List seqMismatchPositionsRcpp(CharacterVector samples,
                           CharacterVector germlines,
                           CharacterVector ignore=CharacterVector::create("N", "-", ".", "?"),
                                bool count_trailing=false) {
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

        positions[i] = wrap(findMismatchPositions(sample, germline, ignore_chars, count_trailing));
    }

    return positions;
}
