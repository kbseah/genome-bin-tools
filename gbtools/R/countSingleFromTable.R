#' Tabulate objects and count how many singletons
#'
#' @param x Object of class data.frame or vector
#' @return Numeric vector of length 1
#' @keywords internal
countSingleFromTable <- function(x) {
    x.tab <- table(x)
    uniq <- length(which(x.tab==1))
    return(uniq)
}
