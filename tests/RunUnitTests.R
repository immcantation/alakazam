#' Lineage tools unit tests
#' @author  Jason Anthony Vander Heiden  
#' @date    2013.11.21

library('RUnit')
PACKAGE.PATH <- '/home/jason/workspace/igpipeline/rigtools_v0.1'
    
source(paste(PACKAGE.PATH, 'lineage/AutomatedPhylipParsimony_v0.1.3.R', sep='/'))
source(paste(PACKAGE.PATH, 'lineage/TreeTopology_v0.1.0.R', sep='/'))

test.suite <- defineTestSuite('lineage', dirs=file.path(paste(PACKAGE.PATH, 'tests', sep='/')), 
                              testFileRegexp = '^Test_\\w+.R')
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)