#include <Rcpp.h>
#include <iostream>
#include <map>
#include <vector>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

namespace {

void buildDistLookupMaps(const NumericMatrix& dist_mat,
                         std::map<std::string, int>& rows_map,
                         std::map<std::string, int>& cols_map) {
    List dist_mat_dims = dist_mat.attr("dimnames");
    CharacterVector dist_mat_rownames = dist_mat_dims[0];
    CharacterVector dist_mat_colnames = dist_mat_dims[1];

    int num_rows = dist_mat_rownames.size();
    int num_cols = dist_mat_colnames.size();

    rows_map.clear();
    cols_map.clear();

    for (int i = 0; i < num_rows; i++) {
        std::string this_row = as<std::string>(dist_mat_rownames[i]);
        if (this_row.empty()) {
            throw std::range_error("Empty row name in dist_mat.");
        }
        rows_map[this_row] = i;
    }

    for (int i = 0; i < num_cols; i++) {
        std::string this_col = as<std::string>(dist_mat_colnames[i]);
        if (this_col.empty()) {
            throw std::range_error("Empty column name in dist_mat.");
        }
        cols_map[this_col] = i;
    }
}

double seqDistWithMaps(const std::string& seq1,
                       const std::string& seq2,
                       const NumericMatrix& dist_mat,
                       const std::map<std::string, int>& rows_map,
                       const std::map<std::string, int>& cols_map) {
    int len_seq1 = seq1.size();
    int len_seq2 = seq2.size();

    if (len_seq1 != len_seq2) {
        throw std::range_error("Sequences of different length.");
    }

    int d_seen = 0;
    int indels = 0;
    double d_sum = 0;

    for (int i = 0; i < len_seq1; i++) {
        int r_idx;
        std::string row_string;
        row_string += static_cast<char>(seq1[i]);
        auto search_row = rows_map.find(row_string);
        if (search_row != rows_map.end()) {
            r_idx = search_row->second;
        } else {
            throw std::range_error("Character not found in dist_mat.");
        }

        int c_idx;
        std::string col_string;
        col_string += static_cast<char>(seq2[i]);
        auto search_col = cols_map.find(col_string);
        if (search_col != cols_map.end()) {
            c_idx = search_col->second;
        } else {
            throw std::range_error("Character not found in dist_mat.");
        }

        double d_i = dist_mat(r_idx, c_idx);

        if (d_i > 0) {
            d_sum = d_sum + d_i;
        } else if ((d_i == -1) & (d_seen != -1)) {
            indels++;
        }
        d_seen = d_i;
    }

    return d_sum + indels;
}

} // namespace


//' Test DNA sequences for equality.
//' 
//' \code{seqEqual} checks if two DNA sequences are identical.
//'
//' @param    seq1    character string containing a DNA sequence.
//' @param    seq2    character string containing a DNA sequence.
//' @param    ignore  vector of characters to ignore when testing for equality.
//'                   Default is to ignore c("N",".","-","?")
//' 
//' @return   Returns \code{TRUE} if sequences are equal and \code{FALSE} if they are not.
//'           Sequences of unequal length will always return \code{FALSE} regardless of
//'           their character values.
//' 
//' @seealso  Used by \link{pairwiseEqual} within \link{collapseDuplicates}.
//'           See \link{seqDist} for calculation Hamming distances between sequences.
//' 
//' @examples
//' # Ignore gaps
//' seqEqual("ATG-C", "AT--C")
//' seqEqual("ATGGC", "ATGGN")
//' seqEqual("AT--T", "ATGGC")
//' 
//' # Ignore only Ns
//' seqEqual("ATG-C", "AT--C", ignore="N")
//' seqEqual("ATGGC", "ATGGN", ignore="N")
//' seqEqual("AT--T", "ATGGC", ignore="N")
//' 
//' @export
// [[Rcpp::export]]
bool seqEqual(std::string seq1, std::string seq2, 
              CharacterVector ignore=CharacterVector::create("N","-",".","?")) {
    
    int ig_len = ignore.length();
    
    ig_len = ignore.length();
    
    int len_seq1 = seq1.length();
    int len_seq2 = seq2.length();
    
    if (len_seq1 != len_seq2) { 
        return (FALSE);
    } else {
        for(int i = 0; i < len_seq1; i++)
        {
            char seq1_char = (char)seq1[i];
            char seq2_char = (char)seq2[i];
            
            if (seq1_char != seq2_char) {
                
                bool ignore_seq1 = FALSE;
                bool ignore_seq2 = FALSE;
                
                for(int j = 0; j < ig_len; j++) {
                    
                    char ig = *(char*)ignore[j];
                    
                    if (ig == seq1_char) {
                        ignore_seq1 = TRUE;
                    }
                    
                    if (ig == seq2_char) {
                        ignore_seq2 = TRUE;
                    }
                }
                if (!ignore_seq1 & !ignore_seq2) {
                    return FALSE;
                }
            }
        }
        return TRUE;
    }  
}


//' Calculate pairwise equivalence between sequences
//' 
//' \code{pairwiseEqual} determined pairwise equivalence between a pairs in a 
//' set of sequences, excluding ambiguous positions (Ns and gaps).
//'
//' @param    seq  character vector containing a DNA sequences.
//'
//' @return   A logical matrix of equivalence between each entry in \code{seq}. 
//'           Values are \code{TRUE} when sequences are equivalent and \code{FALSE}
//'           when they are not.
//' 
//' @seealso  Uses \link{seqEqual} for testing equivalence between pairs.
//'           See \link{pairwiseDist} for generating a sequence distance matrix.
//'           
//' @examples
//' # Gaps and Ns will match any character
//' seq <- c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C", E="NTGGG")
//' d <- pairwiseEqual(seq)
//' rownames(d) <- colnames(d) <- seq
//' d
//' 
//' @export
// [[Rcpp::export]]
LogicalMatrix pairwiseEqual(StringVector seq) {
    
    // allocate the matrix we will return
    LogicalMatrix rmat(seq.length(), seq.length());
    
    for (int i = 0; i < rmat.nrow(); i++) {
        for (int j = 0; j <= i; j++) {
            
            // check seq equal
            std::string row_seq = as<std::string>(seq[i]);
            std::string col_seq = as<std::string>(seq[j]);
            
            bool is_equal = seqEqual(row_seq, col_seq);
            
            // write to output matrix
            rmat(i,j) = is_equal;
            rmat(j,i) = is_equal;
        }
    }
    
    // Add row and column names
    Rcpp::List dimnames = Rcpp::List::create(seq.attr("names"), 
                                             seq.attr("names"));
    rmat.attr("dimnames") = dimnames;
    
    return rmat;
}


// seqDist
// [[Rcpp::export]]
double seqDistRcpp(std::string seq1, std::string seq2, 
                   NumericMatrix dist_mat) {
    std::map<std::string, int> rows_map;
    std::map<std::string, int> cols_map;
    buildDistLookupMaps(dist_mat, rows_map, cols_map);
    return seqDistWithMaps(seq1, seq2, dist_mat, rows_map, cols_map);
}


// pairwiseDist
// [[Rcpp::export]]
NumericMatrix pairwiseDistRcpp(StringVector seq, NumericMatrix dist_mat) {
    // allocate the matrix we will return
    NumericMatrix rmat(seq.length(), seq.length());

    int n_seq = seq.length();
    std::vector<std::string> seq_cpp(n_seq);
    for (int i = 0; i < n_seq; i++) {
        seq_cpp[i] = as<std::string>(seq[i]);
    }

    std::map<std::string, int> rows_map;
    std::map<std::string, int> cols_map;
    buildDistLookupMaps(dist_mat, rows_map, cols_map);
    
    for (int i = 0; i < rmat.nrow(); i++) {
        for (int j = 0; j < i; j++) {
            const std::string& row_seq = seq_cpp[i];
            const std::string& col_seq = seq_cpp[j];

            double distance = seqDistWithMaps(row_seq, col_seq, dist_mat, rows_map, cols_map);
            
            // write to output matrix
            rmat(i,j) = distance;
            rmat(j,i) = distance;
        }
    }
    
    // Add row and column names
    Rcpp::List dimnames = Rcpp::List::create(seq.attr("names"), 
                                             seq.attr("names"));
    rmat.attr("dimnames") = dimnames;
    return rmat;
}


// nonsquareDist
// [[Rcpp::export]]
NumericMatrix nonsquareDistRcpp(StringVector seq, NumericVector indx, NumericMatrix dist_mat)
{
    // defien variables
    int m, n, i, j;
    // extract the sizes. Note: This should be satisfied (n<=m)
    m = indx.size(); //number of rows
    n = seq.size();  //number of columns

    std::vector<std::string> seq_cpp(n);
    for (i = 0; i < n; i++) {
        seq_cpp[i] = as<std::string>(seq[i]);
    }

    // allocate the main matrix
    NumericMatrix rmat(m,n);
    std::fill(rmat.begin(), rmat.end(), NA_REAL);

    std::map<std::string, int> rows_map;
    std::map<std::string, int> cols_map;
    buildDistLookupMaps(dist_mat, rows_map, cols_map);

    // sort and push indices back by 1 to match c++ indexing
    std::sort(indx.begin(), indx.end());
    indx = indx - 1;
    // find the position of the column ids in the indx vector
    NumericVector pos(n);
    for (j = 0; j < n; j++) {
        pos[j] = std::find(indx.begin(), indx.end(), j) - indx.begin();
    }
    // begin filling rmat
    for (i = 0; i < m; i++) {
        int row_id = indx[i];
        const std::string& row_seq = seq_cpp[row_id];
        for (j = 0; j < n; j++) {
            if (!R_IsNA(rmat(i,j))) continue;
            if (row_id == j) rmat(i,j) = 0;
            else {
                const std::string& col_seq = seq_cpp[j];
                rmat(i,j) = seqDistWithMaps(row_seq, col_seq, dist_mat, rows_map, cols_map);
                if (pos[j] < m) rmat(pos[j],indx[i]) = rmat(i,j);
            }
        }
    }
    // Add row and column names
    StringVector subSeq = seq[indx];
    Rcpp::List dimnames = Rcpp::List::create(subSeq.attr("names"),      //rownames
                                             seq.attr("names"));  //colnames
    rmat.attr("dimnames") = dimnames;
    // return matrix
    return rmat;
}
