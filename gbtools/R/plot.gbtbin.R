#' Plot object of class gbtbin
#'
#' Plot GC-coverage or differential coverage plots from a gbt object
#'
#' The plot method for gbt objects can produce both GC-coverage and
#' differential coverage plots. A gbt object contains scaffold data and
#' annotations for a given metagenome, along with the coverage of each
#' scaffold for one or more samples. To plot GC-coverage plots, the sample
#' to use for the coverage data is specified by giving the number of the
#' sample to the slice parameter. To plot differential-coverage plots, the
#' numbers of the two samples for comparison are given to slice as a vector
#' e.g. c(1,2). The first sample is plotted as the x-axis and the second as
#' the y-axis. Supplying a vector with more than 2 elements, or a non-
#' numeric value will return an error.
#' Genome bins can be interactively chosen from a gbt plot with the
#' \code{\link{choosebin}} function; the same slice argument must be passed
#' to the choosebin function as the plot function, otherwise the results will
#' be meaningless!
#' If you wish to plot by user-supplied custom values (added to gbt object by
#' the userAdd() function), indicate which ones to use for x- or y- axis in
#' the userAxis= parameter. For example, to plot user-custom data set 1 in X-
#' and user-custom data set 2 in Y-axis, specify userAxis=c(1,2). To plot vs.
#' GC or coverage values, supply "gc" or "cov" to userAxis, e.g. GC as X-axis
#' and user-value set 2 as Y-axis: userAxis=c("gc",2). If "cov" specified, it
#' will take the coverage set that is given by the slice= parameter (default
#' is the first set of coverage data).
#' 
#' @inheritParams plot.gbt
#' @return New graphics object and plot, or error message
#' @export
#' @seealso \code{\link{gbt}}
#'
plot.gbtbin <- plot.gbt