#' Choose bin from gbt object by defining polygon in coverage-GC or differential coverage plot
#'
#' Choose genome bin by specifying coordinates of a polygon from
#' GC-coverage or differential coverage plot of a gbt object
#'
#' @param x Object of class gbt
#' @param slice Slice parameter for drawing the polygon
#' @param polygon Polygon that would define the bin
#'
#' @importFrom sp point.in.polygon
#' @return Object of class gbtbin
#' @seealso \code{\link{plot.gbt}}
#' @export
#'
choosebinPolygon <- function(x, ... ) UseMethod ("choosebinPolygon")  # Defines generic for choosebinPolygon function
