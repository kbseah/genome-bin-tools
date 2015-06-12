#' Calculate N50 from contig lengths
#'
#' @param x Numeric vector
#' @return N50 value - min contig length that contains half of total bases
#' @keywords internal
getN50 <- function(x) {
    # Adapted from R-bloggers:
    # http://www.r-bloggers.com/calculating-an-n50-from-velvet-output/
    x.sort <- sort(x,decreasing=TRUE)
    n50 <- x[cumsum(x.sort) >= sum(x.sort)/2][1]
    return(n50)
}