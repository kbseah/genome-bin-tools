#' @export
add.gbtbin <- function(x1,x2) {
## Merge two bins; i.e. take their union
    result <- setOperation(x1=x1,x2=x2,shortlist="all")
    result$call[[length(result$call)+1]] <- match.call()  # Record function call 
    return(result)
}
