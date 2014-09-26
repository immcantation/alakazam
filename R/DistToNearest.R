# Generates distance to nearest neighbor by using S5F mutability model
# 
# @author     Namita Gupta, Gur Yaari, Mohamed Uduman
# @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @date       2014.09.24


# Get S5F distance between two sequences of same length broken down into fivemers
#
# @param   seq1   the first nucleotide sequence
# @param   seq2   the second nucleotide sequence
# @return  distance between two sequences based on S5F model
#' @export
dist_seq_fast<-function(seq1,seq2){  
  #Compute distance only on fivemers that have mutations
  fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
  fivemersWithNonNuc <- ( !is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & !is.na(match(substr(seq2,3,3),c("A","C","G","T"))) )
  seq1 <- seq1[fivemersWithMu & fivemersWithNonNuc]
  seq2 <- seq2[fivemersWithMu & fivemersWithNonNuc] 
  a <- tryCatch({
  if(length(seq1)==1){
    #seq1_to_seq2 <- S5F_Substitution[substr(seq2,3,3),seq1] * S5F_Mutability[seq1]
    seq1_to_seq2 <- Targeting[["Targeting"]][substr(seq2,3,3),seq1]
    #seq2_to_seq1 <- S5F_Substitution[substr(seq1,3,3),seq2] * S5F_Mutability[seq2]
    seq2_to_seq1 <- Targeting[["Targeting"]][substr(seq1,3,3),seq2]
  }else{
    seq1_to_seq2 <- sum( diag(Targeting[["Targeting"]][substr(seq2,3,3),seq1]) )
    seq2_to_seq1 <- sum( diag(Targeting[["Targeting"]][substr(seq1,3,3),seq2]) )
  }
  return( mean(c(seq1_to_seq2, seq2_to_seq1)) )
  },error = function(e){
    return(NA)
  })
}


# Given an array of junction sequences, find the distance to the closest sequence
#
# @param   arrJunctions
# @return  distances to the closest sequence
#' @export
getDistanceToClosest <- function(arrJunctions){ 
  
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
      dist_seq_fast(matSequenceFivemers[,i],matSequenceFivemers[,j])
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
#' HS5F model is the SHM targeting model from Yaari, G., et al. Frontiers in Immunology, 2013.
#' MTri model uses the SHM substitution matrix found in Smith, D., et al. J. Immunol., 1996.
#'
#' @param   db          \code{data.frame} which must have the following columns: V_CALL, J_CALL, JUNCTION_GAP_LENGTH, JUNCTION.
#' @param   genotyped   Logical for whether \code{data.frame} is genotyped; if genotyped is true, 
#'                      \code{data.frame} must have column V_CALL_GENOTYPED.
#' @param   first       if TRUE only the first call the gene assignment is used;
#'                      if FALSE the union of ambiguous gene assignments is used to group all sequences with any of those gene calls.
#' @param   model       SHM targeting model, one of c('HS5F','MTri'); see Details.
#' 
#' @return  \code{data.frame} with DIST_NEAREST column added
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "changeo_demo.tab", package="alakazam")
#' df <- readChangeoDb(file)
#' 
#' # Calculate distance to nearest
#' df <- distToNearest(df, genotyped=TRUE, first=FALSE)
#' hist(df$DIST_NEAREST, breaks=50, xlim=c(0,0.02))
#' 
#' @export
distToNearest <- function(db, genotyped=FALSE, first=TRUE, model='HS5F') {  
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
	
	if(model=="HS5F"){
	  load(system.file("extdata", "HS5F_Targeting.RData", package="shm"), envir=.GlobalEnv)
	} else if(model==""){
	  load(system.file("extdata", "MRS5NF_Targeting.RData", package="shm"), envir=.GlobalEnv)
	} else if(model=="MTri"){
	  load(system.file("extdata", "MTri_Targeting.RData", package="shm"), envir=.GlobalEnv)
	} else{
		stop("Error: Unrecognized targeting model\n")
	}
	# Create new column for distance to nearest neighbor
	db$DIST_NEAREST = rep(NA, nrow(db))
	db[,"ROW_ID"] <- 1:nrow(db)
	
	# cat("Calculating distance to nearest neighbor\n")
	db <- arrange( ddply(db, .(V,J,JUNCTION_GAP_LENGTH), mutate, 
                       DIST_NEAREST=getDistanceToClosest(JUNCTION)), 
                 ROW_ID)
	
	db <- db[,!(names(db) %in% c("V","J","ROW_ID"))]
	return(db)
}