#' Subset a gbt or gbtbin object by marker gene taxonomy
#'
#' @param param Taxonomic level to make subset (default "Class")
#' @param value Value of the taxon to make subset (default
#'               "Gammaproteobacteria")
#' @inheritParams winnow
#' @return Object of class gbtbin
#' @seealso \code{\link{winnow}}, \code{\link{gbt}}
#' @export
#'
winnowMark <- function(x,marksource,param,value,save,file) UseMethod("winnowMark")
