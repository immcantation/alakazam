# Obtain CDR3 (or any other region) properties, can be nt or AA 
# 
# @author     Ruoyi Jiang
# @copyright  Copyright 2015 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.2.0
# @date       2015.06.19

#### Constants ####

# Possible regions
VDJ_REGIONS <- c("CDR1", "CDR2", "CDR3", "FWR1", "FWR2", "FWR3")

# Hydropathy index scores
HYDROPATHY <- c(+1.8, -4.5, -3.5, -3.5, +2.5,
                -3.5, -3.5, -0.4, -3.2, +4.5,
                +3.8, -3.9, +1.9, +2.8, -1.6,
                -0.8, -0.7, -0.9, -1.3, +4.2)
names(HYDROPATHY) <- c("A", "R", "N", "D", "C",
                       "Q", "E", "G", "H", "I",
                       "L", "K", "M", "F", "P",
                       "S", "T", "W", "Y", "V")	


#' Translate nucleotide sequences to amino acids
#' 
#' \code{translateDNA} translates nucleotide sequences to AA using functions 
#' from seqinr.
#' 
#' @param   seq   DNA sequence (a string) to be converted to AAs
#' @param   trim  boolean flag to remove 3 nts from both ends of ntseq
#' 
#' @return  string, translated AA stretch
#' 
#' @seealso  \code{\link[seqinr]{translate}}.
#' 
#' @examples
#' library(alakazam)
#' # Load Change-O file
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' translateDNA(df$JUNCTION[1])
#' translateDNA(df$JUNCTION[1], trim=TRUE)
#' translateDNA("ACTGACTCGA")
#' 
#' @export
translateDNA <- function (seq, trim=FALSE) {
	# Returns object seq with 3 nucleotides extracted from either end
	# e.g. "ACTGACTCGA" -> "GACT" (with "ACT" and "CGA" removed)
	if (trim) { seq <- substr(seq, 4, nchar(seq) - 3) }

    # TODO: can replace c2s and s2c within stringi equivalents
    if(nchar(seq) >= 3) {
	    aa <- seqinr::c2s(seqinr::translate(seqinr::s2c(seq)))
	} else {
	    aa <- NA
	}

    return(aa)
}


#' Count patterns
#' 
#' Counts the number of times a "pattern" occurs in "x", a string
#'
#' @param   x            	a string (usually amino acids)
#' @param   pattern         regular expression to be matched in string
#'  
#' @return  numeric, number of times the regular expression was found
#' 
#' @export
countOccurrences <- function(x, pattern) {
  return(sum(gregexpr(pattern, x)[[1]] > 0))
}
			  

#' Calculate GRAVY score
#' 
#' Calculates GRAVY hydrophobicity score from a given AA sequence using the
#' function \code{hydrophobicity} from the \code{Peptides} package.
#'
#' @param   aa_seq          AA sequence (string) to have GRAVY calculated for
#' @param   aa         	 	AA sequence (character vector) from aa_seq with no *
#' 
#' @return  numeric, GRAVY score
#'
#' @seealso  \code{\link[Peptides]{hydrophobicity}}.
#' @examples 
#' aa_seq <- "CARDRSTPWRRGIASTTVRTSW"
#' gravy(aa_seq)
#' gravy(aa_seq, scale= "Aboderin")
#' @export
gravy <- function(aa_seq, scale= "KyteDoolittle"){
    Peptides::hydrophobicity(aa_seq, scale)
}

#' Calculates amino acid properties of regions
#'
#' Obtain data.frame of Region properties. Note this can be used for any region. 
#'
#' @param   region  character vector of either nucleotide or amino acid sequences.
#' @param   nt		boolean, TRUE if the sequences (or sequence) are DNA
#' @param   trim    if \code{TRUE} remove the first and last amino acids from each
#'                  sequence before calculating properties. If \code{FALSE} do
#'                  not modify input sequences.
#' 
#' @return  data.frame, with columns
#' 
#' @examples 
#' region <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
#'             "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
#'             "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
#' regionProperties(region, nt=TRUE, trim=TRUE, region_name="CDR3")
#' 
#' @export
regionProperties <- function(region, nt=FALSE, trim=FALSE, region_name="REGION") {
    #nt=T
    #trim=T
    #region_name=NULL
    
    if (nt) {
        # Check for sequences that are too short
        not_empty <- which(nchar(region) > 2)
        #message(paste(length(not_empty), "Region sequences found."))
        
        # Translate all regions from nt to aa
        region_aa <- rep("", length(region))
        region_aa[not_empty] <- sapply(region[not_empty], translateDNA, trim=trim)
        #message("Region translated to amino acids.")
    } else {
        region_aa <- region
    }
    
    # Calculate region Lengths
    aa_length <- sapply(region_aa, nchar)

    # GRAVY (Grand Average of Hydropathy) index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    # Produces NA if there is a stop codon.
    aa_gravy <- sapply(region_aa, gravy)
    
    # Count the fraction of aa that are positively charged
    aa_positive <- sapply(region_aa, countOccurrences, "[RK]") / aa_length
    
    # Count fraction of aa that are negatively charged
    aa_negative <- sapply(region_aa, countOccurrences, "[DE]") / aa_length
    
    # TODO: Why are we normalizing aliphatic index by length?
    # Aliphatic index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    n_ala <- sapply(region_aa, countOccurrences, "[A]")
    n_val <- sapply(region_aa, countOccurrences, "[V]")
    n_leu_ile <- sapply(region_aa, countOccurrences, "[LI]")
    a <- 2.9
    b <- 3.9
    aa_aliphatic = (n_ala + a*n_val + b*n_leu_ile) / aa_length

    # Count fraction of aa that are aromatic
    aa_aromatic <- sapply(region_aa, countOccurrences, "[FWHY]") / aa_length
    
    # TODO: can be generalized into a loop and individual counts. probably as an argument "aa_count" or something.
    # Count the fraction of aa that are Arg
    aa_arg <- sapply(region_aa, countOccurrences, "[R]") / aa_length
    
    # Count fraction of aa that are His
    aa_his <- sapply(region_aa, countOccurrences, "[H]") / aa_length

    # Count the fraction of aa that are Lys
    aa_lys <- sapply(region_aa, countOccurrences, "[K]") / aa_length
    
    # Count fraction of aa that are Tyr
    aa_tyr <- sapply(region_aa, countOccurrences, "[Y]") / aa_length
    
    # Return the data.frame with amino acid properties
    out_df <- data.frame(AA_LENGTH=aa_length, GRAVY=aa_gravy, 
                         AA_POSITIVE=aa_positive, AA_NEGATIVE=aa_negative, 
                         ALIPHATIC=aa_aliphatic, AROMATIC=aa_aromatic, 
                         ARGININE=aa_arg, HISTIDINE=aa_his, 
                         LYSINE=aa_lys, TYROSINE=aa_tyr)
    colnames(out_df) <- paste0(region_name, "_", colnames(out_df))
    
    return(out_df)
    
}


#' Extract and translate CDR and FWR regions from VDJ region NT sequence
#' Does not extract CDR3!!!
#' 
#' Uses extractVRegion and translate functions to isolate CDR and FWR regions
#' A specific VDJ region can be selected (e.g. region = "CDR1")
#'
#' @param   db 		ChangeoDb, or ChangeoDb column to extract
#'					Must have/be SEQUENCE_IMGT(def.), GERMLINE_IMGT_D_MASK cols
#' @param   col		Column to extract: SEQUENCE_IMGT(def.), GERMLINE_IMGT_D_MASK
#' @param	region	string or vector of strings in VDJ_REGIONS but not "CDR3"
#' 
#' @return  dataframe, with translated regions
#' 
#' @export
translateNonCDR3 <- function(db, col = "SEQUENCE_IMGT", region = VDJ_REGIONS[VDJ_REGIONS != "CDR3"]) {
	output <- NULL
	
	#If db is a single col, the col must be SEQUENCE_IMGT, GERMLINE_IMGT_D_MASK 
	#It can be some other column, but the processing would be nonsensical
	if(is.vector(db)){
		#Remove periods/IMGT gaps
		input <- gsub("[.]","", extractVRegion(db))
	} 
	# If the col that is to be extracted is not in db then error
	else if(!(col %in% colnames(db))){
		stop("Column to be extracted and translated not found")
	}
	#If the col will not lead to the sensical generation of VDJ regions
	else if(!(col %in% c("SEQUENCE_IMGT", "GERMLINE_IMGT_D_MASK"))){
		stop("Column must be SEQUENCE_IMGT or GERMLINE_IMGT_D_MASK")
	}
	# Otherwise, extract the col and generate an input
	else{
		input <- gsub("[.]","", extractVRegion(db[,col]))
	}

	#If the region(s) specified is a subset of CDR1, CDR2, FWR1, FWR2 or FWR3
	if(all(region %in% VDJ_REGIONS[VDJ_REGIONS != "CDR3"])){
		input <- subset(input, select = region)
	} 
	else {
		stop("Region is not CDR1, CDR2, FWR1, FWR2 or FWR3")
	}

	# Translate each column and cbind together
	for(i in 1:ncol(input)){
		output <- cbind(output, unlist(lapply(input[,i], translateDNA))) 
	}
	colnames(output) <- colnames(input)

	return(output)
}


#' Extract and translate CDR3 from junction NT sequence
#'
#' @param   db    	ChangeoDb or vector of strings with nt JUNCTION
#'					ChangeoDb should have column JUNCTION
#'    				
#' @return  vector of strings, contains translated CDR3s
#' 
#' @export
translateCDR3 <- function(db) {
	#If db is a single col, the col must be JUNCTION
	#It can be some other column, but the processing would be nonsensical
	if(is.vector(db)){
		output <- unlist(lapply(db, translateDNA, trim = T))
	}
	#Otherwise, extract the JUNCTION column from the db 
	else{
		if(!("JUNCTION" %in% colnames(db))){
			stop("No JUNCTION column in db")
		}
		output <- unlist(lapply(db[,"JUNCTION"], translateDNA, trim = T))
	}
	return(output)
}


#' Extracts non-CDR3 region properties
#'
#' @param   db 		ChangeoDb or vector of nt or AA strings to extract
#'					ChangeoDb has SEQUENCE_IMGT(def.), GERMLINE_IMGT_D_MASK cols
#'					region must be specified if vector is given
#' @param   col		Column to extract: SEQUENCE_IMGT(def.), GERMLINE_IMGT_D_MASK
#' @param	region		string or vector of strings in VDJ_REGIONS, not "CDR3"
#'						if db is a vector, region will give the region_names
#' @param	nt			boolean, T if the sequences are nt's
#' @param	translate	boolean, T if the output will contain translated db
#' 						functional only if nt = T
#'     					
#' @return  dataframe, contains non-CDR3 region properties (and translations)
#' 
#' @export
extractNonCDR3Properties <- function(db, col = "SEQUENCE_IMGT", region = VDJ_REGIONS[VDJ_REGIONS != "CDR3"], nt = T, translate = T) {
	output <- NULL
	
	#If the user is attempting to translate a non-nucleotide input
	if(!nt & translate){
		stop("Cannot translate non-nucleotide sequences")
	}
	
	#If the db is in nt form...
	#Note db must be or contain the cols SEQUENCE_IMGT, GERMLINE_IMGT_D_MASK.
	#Otherwise the results will be nonsensical 
	if(nt){
		#Add a column to the output that contains the translated regions
		if(translate){
			output <- translateNonCDR3(db, col = col, region = region)
		}
		
		#If the db is a vector SEQUENCE_IMGT, GERMLINE_IMGT_D_MASK...
		if(is.vector(db)){
			input <- subset(gsub("[.]","", extractVRegion(db)), select = region)
			output <- cbind(output, regionProperties(input, nt = T, region_name = region)) 
		} 
		#If db is a dataframe, extract SEQUENCE_IMGT, GERMLINE_IMGT_D_MASK cols.
		else{
			#Remove IMGT gaps and periods and select regions
			input <- subset(gsub("[.]","", extractVRegion(db[,col])), select = region)
			#Generate regionproperties
			for(i in 1:length(colnames(input))){
				output <- cbind(output, regionProperties(input[,colnames(input)[i]], nt = T, region_name = colnames(input)[i]))
			}
		}
	} 
	#Process from the AA output from translate
	else {
		#A string vector db must be a member of CDR1, CDR2, FWR1, FWR2 or FWR3
		if(is.vector(db)){
			input <- db
			output <- cbind(output, regionProperties(input, region_name = region)) 
		} 
		else {	
			#Select regions
			input <- subset(db, select = region)
			#Generate regionproperties
			for(i in 1:length(colnames(input))){
				output <- cbind(output, regionProperties(input[,colnames(input)[i]], region_name = colnames(input)[i])) 
			}
		}
	}
	
	return(output)
}



#' Extracts CDR3 properties
#'  
#' @param   db 		ChangeoDb or vector of nt or AA strings from CDR3 to extract
#'					ChangeoDb has JUNCTION(def.) cols
#' @param   col		Column to extract: SEQUENCE_IMGT(def.), GERMLINE_IMGT_D_MASK
#' @param	region		string or vector of strings in VDJ_REGIONS, not "CDR3"
#' @param	nt			boolean, T if the sequences are nt's
#' @param	translate	boolean, T if the output will contain translated db
#' 						functional only if nt = T
#'     					
#' @return  dataframe, contains CDR3 region properties (and translations)
#' 
#' @export
extractCDR3Properties <- function(db, nt = T, translate = T){
	output <- NULL
	
	if(!nt & translate){
		stop("Cannot translate non-nucleotide sequences")
	}
	
	#If the input is in nt form...
	#Note db must be or contain the cols JUNCTION
	#Otherwise the results will be nonsensical 
	if(nt){
		#Add a column to the output that contains the translated regions
		if(translate){
			output <- cbind(output,CDR3 = translateCDR3(db))
		}
		
		#If the db is a vector of JUNCTION nt, process directly
		if(is.vector(db)){
			input <- unlist(lapply(db, function(x) {substr(x, 4, nchar(x)-3)}))
		}
		#If the db is a dataframe, extract JUNCTION nt and process
		else{
			input <- unlist(lapply(db[,"JUNCTION"], function(x) {substr(x, 4, nchar(x)-3)}))
		}
		
		#Extract regionproperties from the input
		output <- cbind(output, regionProperties(input, nt = T, region_name = "CDR3"))
	}
	#If the db has cols of AA to extract
	else {
		#If the db is a single vector of AA's...
		if(is.vector(db)){
			input <- db
		}
		#If the db is a dataframe, extract the CDR3 column
		else{
			input <- db[,"CDR3"]
		}
		output <- cbind(output,regionProperties(input, region_name = "CDR3"))
	}
	
	return(output)
}