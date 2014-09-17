#' TreeTopology unit tests
#' @author  Jason Anthony Vander Heiden  
#' @date    2013.11.21
#' 
test.stouffer.test <- function()
{
    checkEquals(c(Z=0, p.value=0.5), stouffer.test(c(0, 0)))
    #checkEqualsNumeric(6, factorial(3))
    #checkIdentical(6, factorial(3))
    #checkTrue(2 + 2 == 4, 'Arithmetic works')
    #checkException(log('a'), 'Unable to take the log() of a string')
}

test.deactivation <- function()
{
    DEACTIVATED('Deactivating this test function')
}