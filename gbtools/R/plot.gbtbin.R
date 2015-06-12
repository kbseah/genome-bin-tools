#' Plot object of class gbtbin
#'
#' Plot GC-coverage or differential coverage plots from a gbtbin object
#'
#' See documentation on \code{\link{plot.gbt}} for more details on plotting
#' @inheritParams plot.gbt
#' @seealso \code{\link{plot.gbt}}, \code{\link{points.gbtbin}}
#' @return New graphics window and a plot
#' @export
#'
plot.gbtbin <- function(x, slice="default", ...) {
    ## Defaults to same slice used to choose the bin ################################
    if (slice == "default") { 
        slice <- x$slice
    }
    ## Catch invalid slice values ###################################################
    if (is.na(x$slice) || !is.numeric(x$slice) || length(x$slice) > 2) { 
        cat ("gbtools ERROR: Please specify valid value for slice option for this bin\n")
    }
    ## Else inherit same plot method as gbt class, for simplicty's sake #############
    else {
        plot.gbt (x=x, slice=slice, ...) 
    }
}