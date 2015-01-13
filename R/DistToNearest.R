# Generates distance to nearest neighbor by using S5F mutability model
# 
# @author     Namita Gupta, Gur Yaari, Mohamed Uduman
# @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @date       2014.11.24


# Load targeting model data from the shm package
#
# @param   model   string defining name of the model to load.
#                  One of "hs5f" or "m3n".
# @return  a list containing the substitution (subs) and mutability (mut) models
#' @export
loadModel <- function(model) {
    if (model == "hs5f") { data_file = "HS5F_Targeting.RData" }
    else if (model == "m3n") { data_file = "MTri_Targeting.RData" }
    else { stop("Are you sure you know what you're doing?\n") }
    
    tmp_env <- new.env()
    load(system.file("extdata", data_file, package="shm"), envir=tmp_env)
    Targeting <- get("Targeting", envir=tmp_env)
    rm(tmp_env)
    
    return(list(subs=Targeting[["Substitution"]], mut=Targeting[["Mutability"]]))
}


# Get distance between two sequences of same length broken down into 5-mers.
#
# @param   seq1   first nucleotide sequence.
# @param   seq2   second nucleotide sequence.
# @param   subs    substitution model.
# @param   mut    mutability model.
# @return  distance between two sequences.
#' @export
dist_seq_fast <- function(seq1, seq2, subs, mut){
    #Compute distance only on fivemers that have mutations
    fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
    fivemersWithNonNuc <- ( !is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & !is.na(match(substr(seq2,3,3),c("A","C","G","T"))) )
    seq1 <- seq1[fivemersWithMu & fivemersWithNonNuc]
    seq2 <- seq2[fivemersWithMu & fivemersWithNonNuc]
    a <- tryCatch({
        if(length(seq1)==1){
            seq1_to_seq2 <- subs[substr(seq2,3,3),seq1] * mut[seq1]
            seq2_to_seq1 <- subs[substr(seq1,3,3),seq2] * mut[seq2]
        }else{
            seq1_to_seq2 <- sum( diag(subs[substr(seq2,3,3),seq1]) *  mut[seq1] )
            seq2_to_seq1 <- sum( diag(subs[substr(seq1,3,3),seq2]) *  mut[seq2] )
        }
        return( mean(c(seq1_to_seq2, seq2_to_seq1)) )
    },error = function(e){
        return(NA)
    })
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param   arrJunctions   character vector of junction sequences.
# @param   subs            substitution model.
# @param   mut            mutability model.
# @return  A vector of distances to the closest sequence.
#' @export
getDistanceToClosest <- function(arrJunctions, subs, mut){ 
  
  #Initialize array of distances
  arrJunctionsDist <- rep(NA,length(arrJunctions))
  
  #Filter unique junctions
  arrJunctionsUnique <- unique(arrJunctions)
  
  #Map indexes of unique to its non-unique in the original arrJunctions
  indexJunctions <- match(arrJunctions, arrJunctionsUnique)
  
  #Identify junctions with multiple non-unique sequences and set its distances to 0 
  indexJunctionsCounts <- table(indexJunctions)
  indexRepeated <- as.numeric(names(indexJunctionsCounts)[indexJunctionsCounts>1])
  indexRepeated <- indexJunctions%in%indexRepeated
  arrJunctionsDist[ indexRepeated ] <- rep(0,sum(indexRepeated))
  names(arrJunctionsDist) <- arrJunctions
  
  #Compute distances between junctions
  numbOfUniqueJuctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJuctions)
  if(numbOfUniqueJuctions>1){
    arrJunctionsUnique <- toupper(arrJunctionsUnique)
    arrJunctionsUnique <- gsub('.', '-', arrJunctionsUnique, fixed=TRUE)
    arrJunctionsUnique <- as.vector(sapply(arrJunctionsUnique,function(x){paste("NN",x,"NN",sep="")})) 
    matSequenceFivemers <- sapply(arrJunctionsUnique,
                                  function(x){  
                                    lenString <- nchar(x)
                                    fivemersPos <- 3:(lenString-2)
                                    fivemers <-  substr(rep(x,lenString-4),(fivemersPos-2),(fivemersPos+2))
                                    return(fivemers)
                                  }
                                  , simplify="matrix"
    )
    matDistance <-sapply(1:numbOfUniqueJuctions, function(i)c(rep.int(0,i-1),sapply(i:numbOfUniqueJuctions,function(j){
      dist_seq_fast(matSequenceFivemers[,i],matSequenceFivemers[,j], subs=subs, mut=mut)
    })))
    matDistance <- matDistance + t(matDistance)
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJuctions, function(i){ min(matDistance[-i,i]) })    
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }
  
  #Fill the distances for the sequences that are unique
  arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist[indexJunctionsCounts==1]
  return(round(arrJunctionsDist,4))
}


#' Distance to nearest neighbor
#' 
#' Get distance of every sequence to its nearest sequence sharing same V gene, J gene, and junction length.
#' 
#' hs5f model is the SHM targeting model from Yaari, G., et al. Frontiers in Immunology, 2013.
#' m3n model uses the SHM substitution matrix found in Smith, D., et al. J. Immunol., 1996.
#'
#' @param   db          \code{data.frame} which must have the following columns: V_CALL, J_CALL, JUNCTION_GAP_LENGTH.
#' @param   seq         Field with sequences to compare
#' @param   genotyped   Logical for whether \code{data.frame} is genotyped; if genotyped is true, 
#'                      \code{data.frame} must have column V_CALL_GENOTYPED.
#' @param   first       if TRUE only the first call the gene assignment is used;
#'                      if FALSE the union of ambiguous gene assignments is used to group all sequences with any of those gene calls.
#' @param   model       SHM targeting model, one of c('hs5f','m3n'); see Details.
#' 
#' @return  vector of distances of each sequence to its nearest neighbor
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Calculate distance to nearest
#' df$DIST_NEAREST <- distToNearest(df, genotyped=TRUE, first=FALSE)
#' hist(df$DIST_NEAREST, breaks=50, xlim=c(0,10))
#' 
#' @export
distToNearest <- function(db, seq='JUNCTION', genotyped=FALSE, first=TRUE, model='hs5f') {  
	if(!is.data.frame(db))
		stop('Must submit a data frame')
	
	if(genotyped) { 
		v_col <- "V_CALL_GENOTYPED"
	} else {
		v_col <- "V_CALL"
	}
	j_col <- "J_CALL"
	
	# Parse V and J Column to get gene
	# cat("V+J Column parsing\n")
	
	if(first) {
		db$V <- getGene(db[,v_col])
		db$J <- getGene(db[,j_col])
	} else {
	  db$V <- getGene(db[,v_col], first=FALSE)
	  db$J <- getGene(db[,j_col], first=FALSE)
    # Reassign V genes to most general group of genes
    for(ambig in unique(db$V[grepl(',',db$V)]))
      for(g in strsplit(ambig, split=','))
        db$V2[db$V==g] = ambig
    # Reassign J genes to most general group of genes
	  for(ambig in unique(db$J[grepl(',',db$J)]))
	    for(g in strsplit(ambig, split=','))
	      db$J[db$J==g] = ambig
	}
  
	# Load targeting model
	# cat("Loading Targeting Model\n")
  model_data <- loadModel(model)
  
	# Create new column for distance to nearest neighbor
	db$DIST_NEAREST = rep(NA, nrow(db))
	db[,"ROW_ID"] <- 1:nrow(db)
	
	# cat("Calculating distance to nearest neighbor\n")
	db <- arrange( ddply(db, .(V,J,JUNCTION_GAP_LENGTH), here(mutate), 
                       DIST_NEAREST=getDistanceToClosest(eval(parse(text=seq)), subs=model_data[['subs']], mut=model_data[['mut']])), 
                 ROW_ID)
  
	return(db$DIST_NEAREST)
}