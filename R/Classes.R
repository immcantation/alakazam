# Classes, Generics and Methods

#### Generics ####

# @exportMethod print
#setGeneric("print")

# @exportMethod plot
#setGeneric("plot")


#### Diversity classes ####

#' S4 class defining a clonal abundance curve
#'
#' \code{AbundanceCurve} defines clonal abundance values.
#' 
#' @slot  data    data.frame with relative clonal abundance data and confidence intervals,
#'                containing the following columns:
#'                \itemize{
#'                  \item  \code{GROUP}:  group identifier.
#'                  \item  \code{CLONE}:  clone identifier.
#'                  \item  \code{P}:      relative abundance of the clone.
#'                  \item  \code{LOWER}:  lower confidence inverval bound.
#'                  \item  \code{UPPER}:  upper confidence interval bound.
#'                  \item  \code{RANK}:   the rank of the clone abundance.
#'                }
#' @slot  groups  character vector of group values.
#' @slot  n       numeric vector indication the number of sequences in each group.
#' @slot  nboot   number of bootstrap realizations performed.
#' @slot  ci      confidence interval defining the upper and lower bounds 
#'                (a value between 0 and 1).
#' 
#' @name         AbundanceCurve-class
#' @rdname       AbundanceCurve-class
#' @aliases      AbundanceCurve
#' @exportClass  AbundanceCurve

setClassUnion("charOrNULL",members =c("character","NULL"))
setClassUnion("dfOrNULL",members =c("data.frame","NULL"))
setClassUnion("listOrNULL",members =c("list","NULL"))

setClass("DiversityCalculation", 
         slots=c(
             div="data.frame", 
             abund="dfOrNULL", 
             test="listOrNULL",
             div_group="charOrNULL",
			 div_groups="charOrNULL",
             q="numeric",  
             ci="numeric"))

#' S4 class defining a diversity curve 
#'
#' \code{DiversityCurve} defines diversity (\eqn{D}) scores over multiple diversity 
#' orders (\eqn{Q}).
#' 
#' @slot  data      data.frame defining the diversity curve with the following columns:
#'                  \itemize{
#'                    \item  \code{GROUP}:    group label.
#'                    \item  \code{Q}:        diversity order.
#'                    \item  \code{D}:        mean diversity index over all bootstrap 
#'                                            realizations.
#'                    \item  \code{D_SD}:     standard deviation of the diversity index 
#'                                            over all bootstrap realizations.
#'                    \item  \code{D_LOWER}:  diversity lower confidence inverval bound.
#'                    \item  \code{D_UPPER}:  diversity upper confidence interval bound.
#'                    \item  \code{E}:        evenness index calculated as \code{D} 
#'                                            divided by \code{D} at \code{Q=0}.
#'                    \item  \code{E_LOWER}:  evenness lower confidence inverval bound.
#'                    \item  \code{E_UPPER}:  eveness upper confidence interval bound.
#'                  }
#' @slot  groups    character vector of groups retained in the diversity calculation.
#' @slot  n         numeric vector indication the number of sequences sampled from each group.
#' @slot  nboot     number of bootstrap realizations performed.
#' @slot  ci        confidence interval defining the upper and lower bounds 
#'                  (a value between 0 and 1).
#' 
#' @name         DiversityCurve-class
#' @rdname       DiversityCurve-class
#' @aliases      DiversityCurve
#' @exportClass  DiversityCurve
setClassUnion("dcOrNULL",members =c("DiversityCalculation","NULL"))

setClass("DiversityObject", 
         slots=c(
             data="data.frame",
             bootstrap="data.frame", 
             alpha="dcOrNULL", 
             beta="dcOrNULL", 
             rdi="dcOrNULL", 
             group="character",
             groups="character",
             clone="character",
             status="charOrNULL",
			 copy="charOrNULL", 
             uniform="logical",
             ndepth="numeric", 
			 nboot="numeric", 
             min_n="numeric"))


#' S4 class defining diversity significance
#'
#' \code{DiversityTest} defines the signifance of diversity (\eqn{D}) differences at a 
#' fixed diversity order (\eqn{q}).
#' 
#' @slot  tests    data.frame describing the significance test results with columns:
#'                 \itemize{
#'                   \item  \code{TEST}:        string listing the two groups tested.
#'                   \item  \code{DELTA_MEAN}:  mean of the \eqn{D} bootstrap delta 
#'                                              distribution for the test.
#'                   \item  \code{DELTA_SD}:    standard deviation of the \eqn{D} 
#'                                              bootstrap delta distribution for the test.
#'                   \item  \code{PVALUE}:      p-value for the test.
#'                 }
#' @slot  summary  data.frame containing summary statistics for the diversity index 
#'                 bootstrap distributions, at the given value of \eqn{q}, with columns:
#'                 \itemize{
#'                   \item  \code{GROUP}:   the name of the group.
#'                   \item  \code{MEAN}:    mean of the \eqn{D} bootstrap distribution.
#'                   \item  \code{SD}:      standard deviation of the \eqn{D} bootstrap 
#'                                          distribution.
#'                 }
#' @slot  groups   character vector of groups retained in diversity calculation.
#' @slot  q        diversity order tested (\eqn{q}).
#' @slot  n        numeric vector indication the number of sequences sampled from each group.
#' @slot  nboot    number of bootstrap realizations.
#' 
#' @name         DiversityTest-class
#' @rdname       DiversityTest-class
#' @aliases      DiversityTest
#' @exportClass  DiversityTest
setClass("DiversityTest", 
         slots=c(tests="data.frame",
                 summary="data.frame",
                 groups="character", 
                 q="numeric",
                 n="numeric", 
                 nboot="numeric"))


#### Diversity methods ####

# TODO:  currently not used, implemented as lists
# TODO:  plot method for DiversityTest
# TODO:  summary method for DiversityTest
# TODO:  summary method for DiversityCurve

#' @param    x    AbundanceCurve object
#' 
#' @rdname   AbundanceCurve-class
#' @aliases  AbundanceCurve-method
#' @export
# setMethod("print", c(x="DiversityCalculation"), function(x) { print(x@abund) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotDiversityCurve}.
#' 
#' @rdname   AbundanceCurve-class
#' @aliases  AbundanceCurve-method
#' @export
setMethod("plot", c(x="DiversityCalculation", y="missing"),
          function(x, y, ...) { plotAbundanceCurve(x, ...) })

#' @param    x    DiversityCurve object
#' 
#' @rdname   DiversityCurve-class
#' @aliases  DiversityCurve-method
#' @export
setMethod("print", c(x="DiversityCalculation"), function(x) { print(x@div) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotDiversityCurve}.
#' 
#' @rdname   DiversityCurve-class
#' @aliases  DiversityCurve-method
#' @export
setMethod("plot", c(x="DiversityCalculation", y="missing"),
          function(x, y, ...) { plotDiversityCurve(x, ...) })

#' @param    x    DiversityTest object.
#' 
#' @rdname   DiversityTest-class
#' @aliases  DiversityTest-method
#' @export
setMethod("print", c(x="DiversityCalculation"), function(x) { print(x@test$test) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotDiversityTest}.
#' 
#' @rdname   DiversityTest-class
#' @aliases  DiversityTest-method
#' @export
setMethod("plot", c(x="DiversityCalculation", y="missing"),
          function(x, y, ...) { plotDiversityTest(x, ...) })


#### Lineage classes ####

#' S4 class defining a clone
#' 
#' \code{ChangeoClone} defines a common data structure for perform lineage recontruction
#' from Change-O data.
#' 
#' @slot     data      data.frame containing sequences and annotations. Contains the
#'                     columns \code{SEQUENCE_ID} and \code{SEQUENCE}, as well as any additional 
#'                     sequence-specific annotation columns.
#' @slot     clone     string defining the clone identifier.
#' @slot     germline  string containing the germline sequence for the clone.
#' @slot     v_gene    string defining the V segment gene call.
#' @slot     j_gene    string defining the J segment gene call.
#' @slot     junc_len  numeric junction length (nucleotide count).
#' 
#' @seealso  See \link{makeChangeoClone} and \link{buildPhylipLineage} for use.
#'           
#' @name         ChangeoClone-class
#' @rdname       ChangeoClone-class
#' @aliases      ChangeoClone
#' @exportClass  ChangeoClone
setClass("ChangeoClone", 
         slots=c(data="data.frame",
                 clone="character",
                 germline="character", 
                 v_gene="character", 
                 j_gene="character", 
                 junc_len="numeric"))


#### Topology classes ####

#' S4 class defining edge significance
#'
#' \code{MRCATest} defines the significance of enrichment for annotations appearing at
#' the MRCA of the tree.
#' 
#' @slot  tests         data.frame describing the significance test results with columns:
#'                      \itemize{
#'                        \item  \code{ANNOTATION}:  annotation value.
#'                        \item  \code{COUNT}:       observed count of MRCA positions 
#'                                                   with the given annotation.
#'                        \item  \code{EXPECTED}:    expected mean count of MRCA occurance
#'                                                   for the annotation.
#'                        \item  \code{PVALUE}:      one-sided p-value for the hypothesis that 
#'                                                   the observed annotation abundance is greater 
#'                                                   than expected.
#'                      }
#' @slot  permutations  data.frame containing the raw permutation test data with columns:
#'                      \itemize{
#'                        \item  \code{ANNOTATION}:  annotation value.
#'                        \item  \code{COUNT}:       count of MRCA positions with the 
#'                                                   given annotation.
#'                        \item  \code{ITER}:        numerical index define which 
#'                                                   permutation realization each 
#'                                                   observation corresponds to.
#'                      }
#' @slot  nperm         number of permutation realizations.
#' 
#' @name         MRCATest-class
#' @rdname       MRCATest-class
#' @aliases      MRCATest
#' @exportClass  MRCATest
setClass("MRCATest", 
         slots=c(tests="data.frame",
                 permutations="data.frame",
                 nperm="numeric"))


#' S4 class defining edge significance
#'
#' \code{EdgeTest} defines the significance of parent-child annotation enrichment.
#' 
#' @slot  tests         data.frame describing the significance test results with columns:
#'                      \itemize{
#'                        \item  \code{PARENT}:    parent node annotation.
#'                        \item  \code{CHILD}:     child node annotation
#'                        \item  \code{COUNT}:     count of observed edges with the given 
#'                                                 parent-child annotation set.
#'                        \item  \code{EXPECTED}:  mean count of expected edges for the 
#'                                                 given parent-child relationship.
#'                        \item  \code{PVALUE}:    one-sided p-value for the hypothesis that 
#'                                                  the observed edge abundance is greater 
#'                                                  than expected.
#'                      }
#' @slot  permutations  data.frame containing the raw permutation test data with columns:
#'                      \itemize{
#'                        \item  \code{PARENT}:  parent node annotation.
#'                        \item  \code{CHILD}:   child node annotation
#'                        \item  \code{COUNT}:   count of edges with the given parent-child 
#'                                               annotation set.
#'                        \item  \code{ITER}:    numerical index define which permutation
#'                                               realization each observation corresponds 
#'                                               to.
#'                      }
#' @slot  nperm         number of permutation realizations.
#' 
#' @name         EdgeTest-class
#' @rdname       EdgeTest-class
#' @aliases      EdgeTest
#' @exportClass  EdgeTest
setClass("EdgeTest", 
         slots=c(tests="data.frame",
                 permutations="data.frame",
                 nperm="numeric"))


#### Topology methods ####

#' @param    x    MRCATest object.
#' 
#' @rdname   MRCATest-class
#' @aliases  MRCATest-method
#' @export
setMethod("print", c(x="MRCATest"), function(x) { print(x@tests) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotMRCATest}.
#' 
#' @rdname   MRCATest-class
#' @aliases  MRCATest-method
#' @export
setMethod("plot", c(x="MRCATest", y="missing"),
          function(x, y, ...) { plotMRCATest(x, ...) })

#' @param    x    EdgeTest object.
#' 
#' @rdname   EdgeTest-class
#' @aliases  EdgeTest-method
#' @export
setMethod("print", c(x="EdgeTest"), function(x) { print(x@tests) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotEdgeTest}.
#' 
#' @rdname   EdgeTest-class
#' @aliases  EdgeTest-method
#' @export
setMethod("plot", c(x="EdgeTest", y="missing"),
          function(x, y, ...) { plotEdgeTest(x, ...) })