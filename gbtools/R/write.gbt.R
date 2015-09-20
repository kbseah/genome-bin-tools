#' Export contig list from gbt or gbtbin object
#'
#' @param x Object of class gbt or gbtbin
#' @param file Path to file for writing data (default "gbtools_export")
#'
#' @return External file with list of contigs separated by newlines
#' @seealso \code{\link{gbt}}, \code{\link{gbtbin}}
#' @export
#'

write.gbt <- function(x,  # Object of class gbt or gbtbin
                      file="gbtools_export.list" # File to export list of contigs
                      ) {
    write(as.character(x$scaff$ID),file=file)
}