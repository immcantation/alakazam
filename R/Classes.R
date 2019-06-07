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
#' @slot  abund   data.frame with relative clonal abundance data and confidence intervals,
#'                containing the following columns:
#'                \itemize{
#'                  \item  \code{GROUP}:  group identifier.
#'                  \item  \code{CLONE}:  clone identifier.
#'                  \item  \code{P}:      relative abundance of the clone.
#'                  \item  \code{LOWER}:  lower confidence inverval bound.
#'                  \item  \code{UPPER}:  upper confidence interval bound.
#'                  \item  \code{RANK}:   the rank of the clone abundance.
#'                }
#' @slot  bootstrap data.frame of bootstrapped clonal distributions. 
#' @slot  group     string specifying the name of the group column. 
#' @slot  groups    vector specifying the names of unique groups in group column. 
#' @slot  clone     string specifying the name of the clone column. 
#' @slot  nboot     numeric specifying the number of bootstrap iterations to use.  
#' @slot  uniform   TRUE/FALSE specifying whether or not bootstraps were calculated under rarefaction.
#' @slot  max_n		numeric specifying the number of species used for rarefaction if rarefaction is used. 
#' @slot  min_n		numeric specifying the minumim tolerated number of species in calculation.
#' @slot  ci        confidence interval defining the upper and lower bounds 
#'                  (a value between 0 and 1).
#' 
#' @name         AbundanceCurve-class
#' @rdname       AbundanceCurve-class
#' @aliases      AbundanceCurve
#' @exportClass  AbundanceCurve

setClass("AbundanceCurve", 
         slots=c(
			 abund="data.frame",
             bootstrap="data.frame",
             group="character",
             groups="character",
             clone="character",
             uniform="logical",
             max_n="numeric", 
			 ci="numeric",
			 nboot="numeric", 
             min_n="numeric"))

#' S4 class defining a diversity curve 
#'
#' \code{DiversityCurve} defines diversity (\eqn{D}) scores over multiple diversity 
#' orders (\eqn{Q}).
#' 
#' @slot  div      data.frame defining the diversity curve with the following columns:
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
#' @slot  test      list containing information concerning the result of testing diversity:
#'                  \itemize{
#'                    \item  \code{test}:    data.frame of p-values from testing.
#'                    \item  \code{summary}: data.frame of diversity value means and standard deviations for plotting.
#'                  }
#' @slot  div_group     string specifying the name of the group column in diversity calculation.
#' @slot  div_groups    vector specifying the names of unique groups in group column in diversity calculation.
#' @slot  method         string specifying the type of diversity calculated. 
#' @slot  q         vector of diversity hill diversity indices used for computing diversity.  
#' @slot  ci        confidence interval defining the upper and lower bounds 
#'                  (a value between 0 and 1).
#' 
#' @name         DiversityCurve-class
#' @rdname       DiversityCurve-class
#' @aliases      DiversityCurve
#' @exportClass  DiversityCurve

setClassUnion("listOrNULL",members =c("list","NULL"))

setClass("DiversityCurve", 
          slots=c(
             div="data.frame", 
             test="listOrNULL",
 			 method="character",
             div_group="character",
 			 div_groups="character",
             q="numeric",  
             ci="numeric"))
			  
#### Diversity methods ####

# TODO:  plot method for DiversityTest
# TODO:  summary method for DiversityTest
# TODO:  summary method for DiversityCurve

#' @param    x    AbundanceCurve object
#' 
#' @rdname   AbundanceCurve-class
#' @aliases  AbundanceCurve-method
#' @export
setMethod("print", c(x="AbundanceCurve"), function(x) { print("Abundance Curve") })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotDiversityCurve}.
#' 
#' @rdname   AbundanceCurve-class
#' @aliases  AbundanceCurve-method
#' @export
setMethod("plot", c(x="AbundanceCurve", y="missing"),
          function(x, y, ...) { plotAbundanceCurve(x, ...) })

#' @param    x    DiversityCurve object
#' 
#' @rdname   DiversityCurve-class
#' @aliases  DiversityCurve-method
#' @export
setMethod("print", c(x="DiversityCurve"), function(x) { print(x@method) })

#' @param    y    ignored.
#' @param    ...  arguments to pass to \link{plotDiversityCurve}.
#' 
#' @rdname   DiversityCurve-class
#' @aliases  DiversityCurve-method
#' @export
setMethod("plot", c(x="DiversityCurve", y="missing"),
          function(x, y, ...) { plotDiversityCurve(x, ...) })


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