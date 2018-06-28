#' Subset a gbt or gbtbin object by marker gene taxonomy
#'
#' @param x Object of class gbt or gbtbin
#' @param marksource Which marker set to use (default -- use the first one)
#' @param param Taxonomic level to make subset (default "Class")
#' @param value Value of the taxon to make subset (default "Gammaproteobacteria")
#' @param save Logical: save shortlist of contigs to external file? (default: FALSE)
#' @param file File to save shortlist of contigs
#' 
#' @inheritParams winnow
#' @return Object of class gbtbin
#' @seealso \code{\link{winnow}}, \code{\link{gbt}}
#' @export
#'
winnowMark <- function(x,
                       marksource = NA,
                       param = "Class",
                       value = "Gammaproteobacteria",
                       save = FALSE,
                       file) UseMethod("winnowMark")
