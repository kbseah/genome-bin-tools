#' Take difference between two gbtbin objects
#'
#' Takes the reverse complement of two gbtbin objects. Equivalent to setdiff
#' in R, or left-exclusive-join in SQL. Non commutative!
#'
#' Self explanatory...
#'
#' @inheritParams add
#'
#' @seealso \code{\link{add}}
#'
#' @export
lej <- function(x1, x2) UseMethod ("lej")
