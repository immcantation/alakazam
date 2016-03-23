#include <Rcpp.h>
using namespace Rcpp;

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