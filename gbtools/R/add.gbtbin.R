#' Take union of two gbtbin objects
#'
#' Take union of two gbtbin objects. Equivalent to the R union function
#'
#' Self explanatory...
#' 
#' @param x1 Object of class gbtbin
#' @param x2 Object of class gbtbin
#'
#' @return Object of class gbtbin
#'
#' @seealso \code{\link{lej}}
#' @export
add.gbtbin <- function(x1,x2) {
## Merge two bins; i.e. take their union
    result <- setOperation(x1=x1,x2=x2,shortlist="all")
    result$call[[length(result$call)+1]] <- match.call()  # Record function call 
    return(result)
}
