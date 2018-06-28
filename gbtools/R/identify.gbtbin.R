#' Identify points in gbtools plot
#'
#' Click on GC-coverage plot or differential coverage plot to identify contigs.
#' In cluttered plots it may not be very accurate! Will write the contig ID
#' as a label overlay on the plot.
#'
#' @return Object of class gbtbin containing identified contigs
#' @inheritParams identify.gbt
#'
#' @seealso \code{\link{plot.gbt}}
#' @seealso \code{\link{choosebin}}
#' @export

identify.gbtbin <- identify.gbt 
