#include <Rcpp.h>
#include <iostream>
#include <map>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' Test DNA sequences for equality.
//' 
//' \code{rcpp_testSeqEqual} checks if two DNA sequences are identical.
//'
//' @param    seq1    character string containing a DNA sequence.
//' @param    seq2    character string containing a DNA sequence.
//' @param    ignore  vector of characters to ignore when testing for equality.
//' 
//' @return   Returns \code{TRUE} if sequences are equal and \code{FALSE} if they are not.
//'           Sequences of unequal length will always return \code{FALSE} regardless of
//'           their character values.
//' 
//' @seealso  Used by \link{collapseDuplicates}.
//' 
//' @examples
//' # Ignore gaps
//' rcpp_testSeqEqual("ATG-C", "AT--C")
//' rcpp_testSeqEqual("ATGGC", "ATGGN")
//' rcpp_testSeqEqual("AT--T", "ATGGC")
//' 
//' # Ignore only Ns
//' rcpp_testSeqEqual("ATG-C", "AT--C", ignore="N")
//' rcpp_testSeqEqual("ATGGC", "ATGGN", ignore="N")
//' rcpp_testSeqEqual("AT--T", "ATGGC", ignore="N")
//' @export
// [[Rcpp::export]]
bool rcpp_testSeqEqual(
        std::string seq1, 
        std::string seq2, 
        CharacterVector ignore=CharacterVector() ) {
    
    int ig_len = ignore.length();
    
    // Default ignore characters
    if (ig_len == 0) {
        ignore.push_back("N");
        ignore.push_back("-");
        ignore.push_back(".");
        ignore.push_back("?");
    }
    
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

//' @export
// [[Rcpp::export]]
LogicalMatrix getDistanceMatrix (StringVector rownames) {
    
    // allocate the matrix we will return
    LogicalMatrix rmat(rownames.length(), rownames.length());
    
    for (int i = 0; i < rmat.nrow(); i++) {
        for (int j = 0; j <= i; j++) {
            
            // check seq equal
            std::string row_seq = as<std::string>(rownames[i]);
            std::string col_seq = as<std::string>(rownames[j]);
            
            bool is_equal = rcpp_testSeqEqual(row_seq, col_seq);
            
            // write to output matrix
            rmat(i,j) = is_equal;
            rmat(j,i) = is_equal;
        }
    }
    
    return rmat;
}

//' @export
// [[Rcpp::export]]
IntegerVector validChars (std::string seq1, std::string seq2) {
    
    int len_seq1 = seq1.length();
    int len_seq2 = seq2.length();
    
    if (len_seq1 != len_seq2) {
        throw std::range_error("Sequences of different length.");  
    }
    
    IntegerVector valid_idx=IntegerVector();
    
    for (int i = 0; i < len_seq1; i++)
    {
        char seq1_char = (char)seq1[i];
        char seq2_char = (char)seq2[i];
        
        bool valid_seq1 = ( seq1_char != '.') & (seq1_char != '-');
        bool valid_seq2 = ( seq2_char != '.') & (seq2_char != '-');
        
        if (valid_seq1 | valid_seq2) {
            valid_idx.push_back(i);
        }
    } 
    return valid_idx;
}

/*** R
all.equal(validChars("ATC-C.T", "AT--.TT"), c(0,1,2,4,5,6))
*/

//' Calculate distance between two sequences
//' 
//' \code{getSeqDistance} calculates the distance between two DNA sequences.
//'
//' @param    seq1      character string containing a DNA sequence.
//' @param    seq2      character string containing a DNA sequence.
//' @param    dist_mat  Character distance matrix. Defaults to a Hamming distance 
//'                     matrix returned by \link{getDNAMatrix}. If gap 
//'                     characters, \code{c("-", ".")}, are assigned a value of -1 
//'                     in \code{dist_mat} then contiguous gaps of any run length,
//'                     which are not present in both sequences, will be counted as a 
//'                     distance of 1. Meaning, indels of any length will increase
//'                     the sequence distance by 1. Gap values other than -1 will 
//'                     return a distance that does not consider indels as a special case.
//'
//' @return   Numerical distance between \code{seq1} and \code{seq2}.
//' 
//' @seealso  Nucleotide distance matrix may be built with 
//'           \link{getDNAMatrix}. Amino acid distance matrix may be built
//'           with \link{getAAMatrix}.
//'           
//' @examples
//' # Ungapped examples
//' getSeqDistance("ATGGC", "ATGGG")
//' getSeqDistance("ATGGC", "ATG??")
//' 
//' # Gaps will be treated as Ns with a gap=0 distance matrix
//' getSeqDistance("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=0))
//' 
//' # Gaps will be treated as universally non-matching characters with gap=1
//' getSeqDistance("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=1))
//' 
//' # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
//' getSeqDistance("ATGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))
//' 
//' # Gaps of equivalent run lengths are not counted as gaps
//' getSeqDistance("ATG-C", "ATG-C", dist_mat=getDNAMatrix(gap=-1))
//'
//' # Overlapping runs of gap characters are counted as a single gap
//' getSeqDistance("ATG-C", "AT--C", dist_mat=getDNAMatrix(gap=-1))
//' getSeqDistance("A-GGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))
//' getSeqDistance("AT--C", "AT--C", dist_mat=getDNAMatrix(gap=-1))
//' 
//' # Discontiguous runs of gap characters each count as separate gaps
//' getSeqDistance("-TGGC", "AT--C", dist_mat=getDNAMatrix(gap=-1))
//' 
//' @export
// [[Rcpp::export]]
double rcpp_getSeqDistance(std::string seq1, 
                        std::string seq2, 
                        NumericMatrix dist_mat) {
    
    IntegerVector valid_seq1 = validChars(seq1, seq2);
    
    // seq1 and seq2 have same length
    int len_seqs = valid_seq1.length();
    
    List dist_mat_dims = dist_mat.attr("dimnames");
    //print (dist_mat_dims);
    CharacterVector dist_mat_rownames = dist_mat_dims[0];
    CharacterVector dist_mat_colnames = dist_mat_dims[1];
    int num_rows = dist_mat_rownames.size();
    int num_cols = dist_mat_colnames.size();
    
    List row_key_idx;
    List col_key_idx;
    
    std::map<std::string, int> rows_map;
    std::map<std::string, int> cols_map;
    
    for (int i = 0; i < num_rows; i++)
    {
        //const char *this_col = dist_mat_colnames[i].c_str();
        std::string this_row = as<std::string>(dist_mat_rownames[i]);
        rows_map[this_row] = i;
    }  
    
    for (int i = 0; i < num_cols; i++)
    {
        //const char *this_col = dist_mat_colnames[i].c_str();
        std::string this_col = as<std::string>(dist_mat_colnames[i]);
        cols_map[this_col] = i;
    } 
    
    int d_seen = 0;
    int indels = 0;
    // sum(d[d>0])
    double d_sum = 0;
    
    for (int i = 0; i < len_seqs; i++)
    {
        // find row index
        int row_idx;
        char row_char = (char)seq1[i];
        std::string row_string;
        row_string+=row_char;
        auto search_row = rows_map.find(row_string);
        if(search_row != rows_map.end()) {
            row_idx = search_row->second;
        }
        else {
            throw std::range_error("Character not found in dist_mat.");  
        }
        
        // find col index
        int col_idx;
        char col_char = (char)seq2[i];
        std::string col_string;
        col_string+=col_char;
        auto search_col = cols_map.find(col_string);
        if(search_col != cols_map.end()) {
            col_idx = search_col->second;
        }
        else {
            throw std::range_error("Character not found in dist_mat.");  
        }    
        
        // distance for current i
        double d_i = dist_mat(row_idx, col_idx);
        
        if (d_i > 0){
            // Sum distance
            d_sum = d_sum + d_i;
        } 
        else if ( (d_i == -1 ) &  (d_seen != -1) )
        {
            // Count indel
            indels++;
        }  
        d_seen = d_i;
    }
    
    double distance = d_sum + indels;
    return (distance);
}

//' Calculate pairwise distances between sequences
//' 
//' \code{getSeqMatrix} calculates all pairwise distance between a set of sequences.
//'
//' @param    seq       character vector containing a DNA sequences.
//' @param    dist_mat  Character distance matrix. Defaults to a Hamming distance 
//'                     matrix returned by \link{getDNAMatrix}. If gap 
//'                     characters, \code{c("-", ".")}, are assigned a value of -1 
//'                     in \code{dist_mat} then contiguous gaps of any run length,
//'                     which are not present in both sequences, will be counted as a 
//'                     distance of 1. Meaning, indels of any length will increase
//'                     the sequence distance by 1. Gap values other than -1 will 
//'                     return a distance that does not consider indels as a special case.
//'
//' @return   A matrix of numerical distance between each entry in \code{seq}. 
//'           If \code{seq} is a named vector, row and columns names will be added 
//'           accordingly.
//' 
//' @seealso  Uses \link{getSeqDistance} for calculating distances between pairs.
//'           Nucleotide distance matrix may be built with \link{getDNAMatrix}. 
//'           Amino acid distance matrix may be built with \link{getAAMatrix}. 
//'           
//' @examples
//' # Gaps will be treated as Ns with a gap=0 distance matrix
//' getSeqMatrix(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
//'              dist_mat=getDNAMatrix(gap=0))
//' 
//' # Gaps will be treated as universally non-matching characters with gap=1
//' getSeqMatrix(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
//'              dist_mat=getDNAMatrix(gap=1))
//' 
//' # Gaps of any length will be treated as single mismatches with a gap=-1 distance matrix
//' getSeqMatrix(c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C"), 
//'              dist_mat=getDNAMatrix(gap=-1))
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix rcpp_getSeqMatrix (StringVector rownames, NumericMatrix dist_mat) {
    
    // allocate the matrix we will return
    NumericMatrix rmat(rownames.length(), rownames.length());
    
    for (int i = 0; i < rmat.nrow(); i++) {
        for (int j = 0; j < i; j++) {
            
            // check seq equal
            std::string row_seq = as<std::string>(rownames[i]);
            std::string col_seq = as<std::string>(rownames[j]);
            
            double distance = rcpp_getSeqDistance(row_seq, col_seq, dist_mat);
            
            // write to output matrix
            rmat(i,j) = distance;
            rmat(j,i) = distance;
        }
    }
    
    Rcpp::List dimnames = Rcpp::List::create(rownames.attr("names"), 
                                             rownames.attr("names"));
    rmat.attr("dimnames") = dimnames;
    return rmat;
}