#' Subset a gbt or gbtbin object by GC, length, or coverage cutoffs
#'
#' @param x Object of class gbt or gbtbin
#' @param gc Vector giving min and max GC values (default c(0,Inf))
#' @param len Vector giving min and max contig lengths (default c(0,Inf))
#' @param covmin Vector giving minimum coverage values per sample (default NA)
#' @param covmax Vector giving maximum coverage values per sample (default NA)
#' @param slice Vector of sample numbers for the coverage cutoffs (default NA)
#' @param save Save list of contigs to external file? (logical, defautl FALSE)
#' @param file File name for export of contig list.
#' @return Object of class gbtbin
#' @seealso \code{\link{winnowMark}}, \code{\link{gbt}}
#' @export
#'
winnow.gbtbin <- winnow.gbt  # Inherit behavior of gbt method
