#' Perform connectivity fishing with Fastg and paths files from SPAdes 3.6.2+
#'
#' Find contigs connected to a given gbtbin and report a new gbtbin object
#' Requires Fastg and paths files produced by SPAdes 3.6.2+; old version of
#' SPAdes produces different format.
#'
#' @param x Object of class gbt (parent object of the gbtbin object)
#' @param bin Object of class gbtbin, derived from x above
#' @param fastg.file Fastg formatted assembly graph from SPAdes
#' @param paths Paths file mapping assembly graph edge names to scaffold/contig names
#' @param save Logical: Save list of fished contigs to external file? (Default: No)
#' @param file File name to save list of fished contigs, if save=TRUE
#' @return Object of class gbtbin
#' @export

fastgFish <- function (x, bin, fastg.file, paths.file, save, file) UseMethod ("fastgFish") # Define generic for fastgFish function