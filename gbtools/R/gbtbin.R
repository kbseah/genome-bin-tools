#' Create gbtbin from gbt object
#'
#' Creates object of class gbtbin from existing gbt object, representing a
#' genome bin defined from a metagenome assembly.
#'
#' This function is called by \code{\link{choosebin}}, but can be called
#' directly by the user, given a shortlist of contigs in the bin.
#'
#' @param shortlist List of contigs in the new bin (required)
#' @param x Object of class gbt
#' @param slice Numeric vector describing which coverage values were used to
#'               draw the plot from which the bin was selected by choosebin.
#'               Use NA if not calling this function from choosebin.
#' @param taxon For internal use
#' @param points Numeric vector describing points used to select the bin from
#'                a plot by choosebin. Use NA if not calling this function
#'                from choosebin.
#' @param save Save list of contigs to external file? (logical, default FALSE)
#' @param file File name to save contigs, if save=TRUE.
#'
#' @return Object of class gbtbin
#' @seealso \code{\link{gbt}}, \code{\link{choosebin}}
#' @export
gbtbin <- function(shortlist,x,slice,taxon,points,save,file) UseMethod("gbtbin")
