#' Generic operation for merging and subsetting two gbtbin objects
#' 
#' @param x1 Object of class gbtbin
#' @param x2 Object of class gbtbin
#' @param shortlist Vector of contig IDs to make new bin
#' @return Object of class gbtbin
#' @keywords internal

setOperation <- function(x1, x2, shortlist) UseMethod("setOperation")
