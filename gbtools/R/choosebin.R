#' Choose bin interactively from plot of gbt object
#'
#' Choose genome bin from GC-coverage or differential coverage plot of a
#' gbt object
#'
#' @param x Object of class gbt, used to generate the plot
#' @param slice The same slice parameter used to generate the plot
#' @inheritParams pickBinPoints
#'
#' @importFrom sp point.in.polygon
#' @return Object of class gbtbin
#' @seealso \code{\link{plot.gbt}}
#' @export
#'
choosebin <- function(x, ... ) UseMethod ("choosebin")  # Defines generic for choosebin function
